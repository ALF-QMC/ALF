#!/usr/bin/env python3

import numpy as np
import h5py
import time
import logging
from multiprocessing import Pool, current_process
import argparse
import os
import shutil

# ---------- 日志配置 ----------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("hf_scan.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("HF_Scan")

TMP_DIR = "tmp_results"
FINAL_H5 = "tV_HF_results.h5"

# ===== 合并函数 =====
def merge_h5_files(tmp_dir, final_file):
    with h5py.File(final_file, "a") as fout:
        for fname in os.listdir(tmp_dir):
            if fname.endswith(".h5"):
                with h5py.File(os.path.join(tmp_dir, fname), "r") as fin:
                    for key in fin.keys():
                        if key not in fout:
                            fin.copy(key, fout)
    logger.info(f"Merge completed! Final results saved in {final_file}")

    shutil.rmtree(tmp_dir)
    logger.info(f"Temporary directory '{tmp_dir}' removed.")


# ---------- 工具函数 ----------
def get_index(ix, iy, sublattice, Lx, Ly):
    return 2 * (ix * Ly + iy) + sublattice

def build_checkerboard_neighbors(Lx, Ly, t1=1.0, t2p=0.2, pbc=True):
    rows, cols, values = [], [], []
    for ix in range(Lx):
        for iy in range(Ly):
            jb = get_index(ix, iy, 1, Lx, Ly)
            offsets = [(0, 0), (0, 1), (1, 0), (1, 1)]
            for dx, dy in offsets:
                ia_x = (ix + dx) % Lx if pbc else ix + dx
                ia_y = (iy + dy) % Ly if pbc else iy + dy
                if 0 <= ia_x < Lx and 0 <= ia_y < Ly:
                    ia = get_index(ia_x, ia_y, 0, Lx, Ly)
                    rows += [ia, jb]
                    cols += [jb, ia]
                    values += [-t1, -t1]
    for ix in range(Lx):
        for iy in range(Ly):
            for dx, dy, sign in [(1, 0, 1.0), (0, 1, -1.0)]:
                ix2 = (ix + dx) % Lx if pbc else ix + dx
                iy2 = (iy + dy) % Ly if pbc else iy + dy
                if 0 <= ix2 < Lx and 0 <= iy2 < Ly:
                    for s in [0, 1]:
                        i1 = get_index(ix, iy, s, Lx, Ly)
                        i2 = get_index(ix2, iy2, s, Lx, Ly)
                        factor = sign if s == 0 else -sign
                        rows += [i1, i2]
                        cols += [i2, i1]
                        values += [-factor * t2p, -factor * t2p]
    return np.array(rows), np.array(cols), np.array(values)

def build_bond_mapping(Lx, Ly, pbc=True):
    bond_list = []
    for ix in range(Lx):
        for iy in range(Ly):
            jb = get_index(ix, iy, 1, Lx, Ly)
            offsets = [(0, 0), (0, 1), (1, 0), (1, 1)]
            for dx, dy in offsets:
                ia_x = (ix + dx) % Lx if pbc else ix + dx
                ia_y = (iy + dy) % Ly if pbc else iy + dy
                if 0 <= ia_x < Lx and 0 <= ia_y < Ly:
                    ia = get_index(ia_x, ia_y, 0, Lx, Ly)
                    bond_list.append((ia, jb))
    return np.array(bond_list, dtype=np.int32)

def build_bond_mapping_v2(Lx, Ly, pbc=True):
    bond_list = []
    for ix in range(Lx):
        for iy in range(Ly):
            for dx, dy in [(1, 0), (0, 1)]:
                ix2 = (ix + dx) % Lx if pbc else ix + dx
                iy2 = (iy + dy) % Ly if pbc else iy + dy
                if 0 <= ix2 < Lx and 0 <= iy2 < Ly:
                    ia1 = get_index(ix, iy, 0, Lx, Ly)
                    ia2 = get_index(ix2, iy2, 0, Lx, Ly)
                    bond_list.append((ia1, ia2))
                if 0 <= ix2 < Lx and 0 <= iy2 < Ly:
                    ib1 = get_index(ix, iy, 1, Lx, Ly)
                    ib2 = get_index(ix2, iy2, 1, Lx, Ly)
                    bond_list.append((ib1, ib2))
    return np.array(bond_list, dtype=np.int32)

#def build_hamiltonian(n_expect, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite):
#    H = np.zeros((Nsite, Nsite))
#    H[rows, cols] = values
#    for ia, ib in bond_pairs1:
#        H[ia, ia] += V1 * n_expect[ib] - V1 * 0.5
#        H[ib, ib] += V1 * n_expect[ia] - V1 * 0.5
#    for i1, i2 in bond_pairs2:
#        H[i1, i1] += V2 * n_expect[i2] - V2 * 0.5
#        H[i2, i2] += V2 * n_expect[i1] - V2 * 0.5
#    return H

def build_hamiltonian(n_expect, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite, Lx, Ly, h_pin=1e-5):
    H = np.zeros((Nsite, Nsite))
    H[rows, cols] = values
    for ia, ib in bond_pairs1:
        H[ia, ia] += V1 * n_expect[ib] - V1 * 0.5
        H[ib, ib] += V1 * n_expect[ia] - V1 * 0.5
    for i1, i2 in bond_pairs2:
        H[i1, i1] += V2 * n_expect[i2] - V2 * 0.5
        H[i2, i2] += V2 * n_expect[i1] - V2 * 0.5

    # -------- Pinning Field on boundary --------
    for iy in range(Ly):
        for sublattice in [0, 1]:   # 0: a, 1: b
            idx = get_index(0, iy, sublattice, Lx, Ly)
            sign = -1 if sublattice == 0 else 1
            H[idx, idx] += h_pin * sign
    return H

def total_density_matrix(H, Nelec):
    eigvals, eigvecs = np.linalg.eigh(H)
    F_occ = np.zeros(H.shape[0])
    F_occ[:Nelec] = 1.0
    n_matrix = eigvecs @ np.diag(F_occ) @ eigvecs.T.conj()
    slater_matrix = eigvecs[:, :Nelec]   # N_site x N_elec
    return n_matrix, slater_matrix

def get_random_key(seed=None):
    return np.random.default_rng(seed if seed is not None else int(time.time()) % (2**32 - 1))

def generate_initial_state(state_type, Lx, Ly, rng, amp_params=(0.3,0.7), noise_level=0.08):
    Nsite = 2 * Lx * Ly
    n_init = np.zeros(Nsite)
    if state_type == "y-stripe":
        for ix in range(Lx):
            for iy in range(Ly):
                val = amp_params[1] if (iy % 2 == 0) else amp_params[0]
                for s in [0, 1]:
                    idx = get_index(ix, iy, s, Lx, Ly)
                    n_init[idx] = val
    elif state_type == "diag-stripe":
        for ix in range(Lx):
            for iy in range(Ly):
                val = amp_params[1] if ((ix + iy) % 2 == 0) else amp_params[0]
                for s in [0, 1]:
                    idx = get_index(ix, iy, s, Lx, Ly)
                    n_init[idx] = val
    elif state_type == "diag-stripe-complex":
        high_main, low_main = amp_params  # 例如 (0.7, 0.3)
        for ix in range(Lx):
            for iy in range(Ly):
                # 第一层：决定当前点属于高条纹还是低条纹
                if (ix + iy) % 2 == 0:
                    base_val = high_main
                else:
                    base_val = low_main
                # 第二层：在对角线上交错强弱
                offset = 0.05  # 交错幅度差
                if (ix % 2) == 0:
                    val = base_val + offset
                else:
                    val = base_val - offset
                for s in [0, 1]:
                    idx = get_index(ix, iy, s, Lx, Ly)
                    n_init[idx] = val
    elif state_type == "ab-cdw":
        for ix in range(Lx):
            for iy in range(Ly):
                for s in [0, 1]:
                    val = amp_params[1] if s == 0 else amp_params[0]
                    idx = get_index(ix, iy, s, Lx, Ly)
                    n_init[idx] = val
    elif state_type == "random":
        n_init = rng.uniform(0.0, 1.0, size=Nsite)
    else:
        raise ValueError(f"Unknown initial state type: {state_type}")
    noise = noise_level * rng.normal(size=Nsite)
    n_init = np.clip(n_init + noise, 0.0, 1.0)
    return n_init

def run_optimization_with_ninit(n_init, Lx, Ly, t1, t2p, V1, V2, Nelec, maxiter=500, pbc=True, alpha=0.5):
    Nsite = int(2 * Lx * Ly)
    bond_pairs1 = build_bond_mapping(Lx, Ly, pbc)
    bond_pairs2 = build_bond_mapping_v2(Lx, Ly, pbc)
    rows, cols, values = build_checkerboard_neighbors(Lx, Ly, t1, t2p, pbc)

    n_expect = n_init
    for i in range(maxiter):
        pin_idx = get_index(0, 0, 0, Lx, Ly)
        #H = build_hamiltonian(n_expect, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite)
        H = build_hamiltonian(n_expect, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite, Lx, Ly, h_pin=1e-5)
        n_matrix, slater_matrix = total_density_matrix(H, Nelec)
        flat_n_new = np.real(np.diag(n_matrix))
        n_expect_next = (1 - alpha) * n_expect + alpha * flat_n_new
        diff = np.linalg.norm(n_expect_next - n_expect)
        if diff < 1e-7:
            break
        n_expect = n_expect_next

    n_final = n_expect.reshape((Lx, Ly, 2))

    interaction1 = V1 * np.sum(n_expect[bond_pairs1[:,0]] * n_expect[bond_pairs1[:,1]])
    interaction2 = V2 * np.sum(n_expect[bond_pairs2[:,0]] * n_expect[bond_pairs2[:,1]])
    const_term = (V1/4) * bond_pairs1.shape[0] + (V2/4) * bond_pairs2.shape[0]
    energy = np.trace(H @ n_matrix).real - 0.5 * (interaction1 + interaction2) + const_term

    return n_final, energy, n_expect, slater_matrix 

def compute_staggered_CDW(n_mean):
    Lx, Ly, _ = n_mean.shape
    total = 0.0
    for ix in range(Lx):
        for iy in range(Ly):
            for s in [0, 1]:
                #idx_factor = (-1) ** (ix + iy + s)
                idx_factor = (-1) ** (s)
                total += idx_factor * (n_mean[ix, iy, s] - 0.5)
    Nsite = Lx * Ly * 2
    return total / Nsite

# ===== 核心计算函数：检测最终文件并写独立临时文件 =====
def worker_task(args):
    v1, v2, params_base, initial_state_configs, random_keys, Lx, Ly = args
    group_key = f"V1_{v1:.2f}_V2_{v2:.2f}"

    proc_name = current_process().name

    # 优先检测最终文件是否已有结果
    if os.path.exists(FINAL_H5):
        with h5py.File(FINAL_H5, "r") as final_h5:
            if group_key in final_h5:
                logger.info(f"{proc_name} | {group_key} already in final file. Skipping.")
                return

    # 检查并创建临时文件
    tmp_file = os.path.join(TMP_DIR, f"{proc_name}.h5")
    if not os.path.exists(tmp_file):
        with h5py.File(tmp_file, "w"): pass

    with h5py.File(tmp_file, "a") as h5file:
        if group_key in h5file:
            logger.info(f"{proc_name} | {group_key} exists in tmp. Skipping.")
            return

    # ===== 计算过程 =====
    best_energy = 1e10
    best_result = {}

    for state_type, amp_list in initial_state_configs:
        if state_type == "random":
            keys = random_keys
            amps = [None] * len(random_keys)
        else:
            keys = [get_random_key() for _ in amp_list]
            amps = amp_list

        for key, amp in zip(keys, amps):
            n_init = generate_initial_state(state_type, Lx, Ly, key, amp_params=amp if amp else (0.3,0.7))
            n_mean, energy, n_expect, slater = run_optimization_with_ninit(
                n_init, Lx, Ly, t1=params_base['t1'], t2p=params_base['t2p'],
                V1=v1, V2=v2, Nelec=params_base['Nelec'], maxiter=params_base['maxiter'], pbc=params_base['pbc']
            )
            if energy < best_energy:
                best_energy = energy
                best_result = dict(n_mean=n_mean, slater=slater)

    stagger_cdw = compute_staggered_CDW(best_result['n_mean'])

    with h5py.File(tmp_file, "a") as h5file:
        grp = h5file.create_group(group_key)
        grp.create_dataset("n_mean", data=best_result['n_mean'])
        grp.create_dataset("energy", data=best_energy)
        grp.create_dataset("stagger_cdw", data=stagger_cdw)
        grp.create_dataset("slater", data=best_result['slater'])

    logger.info(f"{proc_name} | Saved: {group_key} | Energy = {best_energy:.6f} | CDW = {stagger_cdw:.6f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nproc", type=int, default=2, help="Number of processes")
    args = parser.parse_args()

    Lx, Ly = 4, 4
    params_base = dict(t1=1.0, t2p=0.5, Nelec=16, maxiter=500, pbc=True)

    V1_list = np.linspace(0.1, 2, 10)
    V2_list = [0.00]

    initial_state_configs = [
        ("y-stripe",  [(0.3, 0.7),(0.4,0.6)]),
        ("diag-stripe", [(0.3, 0.7), (0.4,0.6)]),
        ("diag-stripe-complex", [(0.7, 0.3), (0.4, 0.6)]),
        ("ab-cdw",   [(0.2, 0.8), (0.4,0.6)]),
        ("random",   [])
    ]

    num_random_trials = 5
    random_keys = [get_random_key(seed=100 + i*32148975167) for i in range(num_random_trials)]

    os.makedirs(TMP_DIR, exist_ok=True)

    tasks = []
    for v1 in V1_list:
        for v2 in V2_list:
            tasks.append( (v1, v2, params_base, initial_state_configs, random_keys, Lx, Ly) )

    logger.info(f"Total tasks: {len(tasks)} | Using {args.nproc} processes.")

    with Pool(processes=args.nproc) as pool:
        pool.map(worker_task, tasks)

    merge_h5_files(TMP_DIR, FINAL_H5)

    logger.info("All tasks completed!")
