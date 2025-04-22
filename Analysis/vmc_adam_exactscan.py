#!/usr/bin/env python3
import jax
import jax.numpy as jnp
import numpy as np
import optax
import h5py
from jax.scipy.linalg import eigh
from joblib import Parallel, delayed
import os
import csv
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

jax.config.update("jax_enable_x64", True)

###concurrent.futures
#mp.set_start_method('spawn', force=True)
#def parallel_run(batch, batch_size):
#    with ProcessPoolExecutor(max_workers=batch_size) as executor:
#        results = list(executor.map(vmc_trial_worker, batch))
#    return results

##multiprocessing
def parallel_run(batch, batch_size):
    ctx = mp.get_context('spawn')
    with ctx.Pool(processes=batch_size) as pool:
        results = pool.map(vmc_trial_worker, batch)
    return results

# -----------------------------------------
# Utils
# -----------------------------------------
def pack(M):
    return jnp.concatenate([M.real.flatten(), M.imag.flatten()])

def unpack(p, N, Nf):
    M_real = p[:N*Nf].reshape((N, Nf))
    M_imag = p[N*Nf:].reshape((N, Nf))
    return M_real + 1j * M_imag

def compute_G(M):
    Q, _ = jnp.linalg.qr(M)
    return Q @ Q.conj().T

# Safe HDF5 save function
def save_phi_to_h5_safe(filename, key, M):
    M = np.asarray(M, dtype=np.complex128)
    data = np.stack([M.real, M.imag], axis=-1)
    with h5py.File(filename, "a") as f_h5:
        if key in f_h5:
            print(f"HDF5 key '{key}' already exists, skipping write.")
        else:
            f_h5.create_dataset(key, data=data, dtype=np.float64)

# write CSV
def append_to_csv(filename, header, row):
    file_exists = os.path.exists(filename)
    mode = 'a' if file_exists else 'w'
    with open(filename, mode, newline="") as f_csv:
        writer = csv.writer(f_csv)
        if not file_exists:
            writer.writerow(header)
        writer.writerow(row)

def save_partial_result(loss, params, label):
    fname = f"trial_{label.replace('=','').replace('.','p')}.npz"
    np.savez(fname, loss=loss, params=params)

def compute_staggered_order(G, sublattices):
    signs = jnp.where(sublattices == 0, 1.0, -1.0)
    G_diag = jnp.real(jnp.diag(G))
    return jnp.sum(signs * G_diag)

# -----------------------------------------
# Lattice Construction
# -----------------------------------------
def make_checkerboard_lattice(Lx, Ly):
    N_sites = 2 * Lx * Ly
    coords = []
    sublattices = []
    for y in range(Ly):
        for x in range(Lx):
            coords.append((x, y))
            sublattices.append(0)
            coords.append((x + 0.5, y + 0.5))
            sublattices.append(1)
    return N_sites, jnp.array(coords), jnp.array(sublattices)

def get_bond_lists(Lx, Ly):
    def site_index(x, y, sub):
        return 2 * ((y % Ly) * Lx + (x % Lx)) + sub

    nn_bonds, nn_dirs, nnn_bonds, nnn_dirs = [], [], [], []

    for y in range(Ly):
        for x in range(Lx):
            ia = site_index(x, y, 0)
            ib = site_index(x, y, 1)
            nn_bonds += [(ib, ia), (ib, site_index(x, y+1, 0)),
                         (ib, site_index(x+1, y, 0)), (ib, site_index(x+1, y+1, 0))]
            nn_dirs += [(1, 'n1'), (-1, 'n1'), (-1, 'n1'), (1, 'n1')]

            ja_x, ja_y = site_index(x+1, y, 0), site_index(x, y+1, 0)
            jb_x, jb_y = site_index(x+1, y, 1), site_index(x, y+1, 1)
            nnn_bonds += [(ia, ja_x), (ia, ja_y), (ib, jb_x), (ib, jb_y)]
            nnn_dirs += [(0, 'x'), (0, 'y'), (1, 'x'), (1, 'y')]

    return list(set(nn_bonds)), nnn_bonds, nn_dirs, nnn_dirs

# -----------------------------------------
# Hamiltonian and Energy
# -----------------------------------------
def make_free_hamiltonian(N, nn_bonds, nnn_bonds, nn_dirs, nnn_dirs, sublattices, t1, t2, m):
    H = jnp.zeros((N, N), dtype=jnp.complex128)
    for (i, j), (sub, _) in zip(nn_bonds, nn_dirs):
        phase = -t1 * jnp.exp(sub * jnp.pi / 4.0)
        H = H.at[i, j].set(phase)
        H = H.at[j, i].set(jnp.conj(phase))
    for (i, j), (sub, direction) in zip(nnn_bonds, nnn_dirs):
        sign = -t2 if (sub == 0 and direction == 'x') or (sub == 1 and direction == 'y') else t2
        H = H.at[i, j].set(sign)
        H = H.at[j, i].set(sign)
    for i in range(N):
        H = H.at[i, i].add(-m if sublattices[i] == 0 else m)
    return H

def energy_wick(G, nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2):
    E_kin, E_int = 0.0, 0.0
    for i, j in nn_bonds:
        E_kin += -t1 * (G[i, j] + G[j, i])
    for (i, j), (sub, direction) in zip(nnn_bonds, nnn_dirs):
        sign = -t2 if (sub == 0 and direction == 'x') or (sub == 1 and direction == 'y') else t2
        E_kin += sign * (G[i, j] + G[j, i])
    for (i, j), V in zip(nn_bonds + nnn_bonds, [V1]*len(nn_bonds) + [V2]*len(nnn_bonds)):
        ni, nj = G[i, i], G[j, j]
        E_int += V * ((ni - 0.5)*(nj - 0.5) - G[i, j]*G[j, i])
    return jnp.real(E_kin + E_int)

# -----------------------------------------
# Worker for parallel
# -----------------------------------------
def vmc_trial_worker(args):
    label, M_init, N, Nf, nn_bonds, nnn_bonds, nnn_dirs, t1, t2, V1, V2, lr, n_iter = args
    params = pack(M_init)
    opt = optax.adam(lr)
    opt_state = opt.init(params)

    @jax.jit
    def step(params, opt_state):
        M = unpack(params, N, Nf)
        G = compute_G(M)
        loss = energy_wick(G, nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2)
        grads = jax.grad(lambda p: energy_wick(compute_G(unpack(p, N, Nf)), nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2))(params)
        updates, opt_state = opt.update(grads, opt_state)
        return optax.apply_updates(params, updates), opt_state, loss

    best_loss = jnp.inf
    best_params = None
    for _ in range(n_iter):
        params, opt_state, loss = step(params, opt_state)
        if loss < best_loss:
            best_loss = loss
            best_params = params

    #save_partial_result(best_loss, best_params, label)
    return float(best_loss), np.array(best_params), label

def run_trials_with_refinement(
    Lx, Ly, t1, t2, V1, V2, rough_iter, refine_iter, top_k, batch_size, lr
):
    N, coords, sublattices = make_checkerboard_lattice(Lx, Ly)
    Nf = N // 2
    nn_bonds, nnn_bonds, nn_dirs, nnn_dirs = get_bond_lists(Lx, Ly)
    trials = []

    m_list = jnp.linspace(0.01, 1.01, 21)
    seed_list = range(0, 260, 13)

    for m in m_list:
        H0 = make_free_hamiltonian(N, nn_bonds, nnn_bonds, nn_dirs, nnn_dirs, sublattices, t1, t2, m)
        evals, evecs = eigh(H0)
        M_init = evecs[:, :Nf]
        trials.append(("m=%.3f" % m, M_init))

    for seed in seed_list:
        np.random.seed(seed)
        M_real = np.random.randn(N, Nf)
        M_imag = np.random.randn(N, Nf)
        M_init = M_real + 1j * M_imag
        trials.append(("seed=%d" % seed, M_init))

    task_args = [
        (label, M_init, N, Nf, nn_bonds, nnn_bonds, nnn_dirs, t1, t2, V1, V2, lr, rough_iter)
        for label, M_init in trials
    ]

    rough_results = []
    for i in range(0, len(task_args), batch_size):
        batch = task_args[i:i+batch_size]
        #results = Parallel(n_jobs=batch_size)(delayed(vmc_trial_worker)(args) for args in batch)
        results = parallel_run(batch, batch_size)
        rough_results.extend(results)

    sorted_results = sorted(rough_results, key=lambda x: x[0])
    top_results = sorted_results[:top_k]

    refined = []
    for loss, params_arr, label in top_results:
        M_init = unpack(jnp.array(params_arr), N, Nf)
        refined_loss, refined_params, _ = vmc_trial_worker(
            (label + "_refined", M_init, N, Nf, nn_bonds, nnn_bonds, nnn_dirs, t1, t2, V1, V2, lr, refine_iter)
        )
        refined.append((refined_loss, refined_params, label + "_refined"))

    best_loss, best_params, best_label = min(refined, key=lambda x: x[0])
    M_best = unpack(jnp.array(best_params), N, Nf)
    return M_best

# -----------------------------------------
def sweep_v1(t1, t2, V2, V1_list, mode="append", **kwargs):
    Lx, Ly = kwargs['Lx'], kwargs['Ly']
    N, coords, sublattices = make_checkerboard_lattice(Lx, Ly)

    csv_file = "results_summary.csv"
    h5_file = "all_phi_trials.h5"
    header = ["t1", "t2", "V1", "V2", "energy", "staggered_order", "filename"]

    # 读取已完成的 V1 列表，避免重复计算
    done_V1 = set()
    if os.path.exists(csv_file) and mode == "append":
        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if float(row["t1"]) == t1 and float(row["t2"]) == t2 and float(row["V2"]) == V2:
                    done_V1.add(float(row["V1"]))

    for V1 in V1_list:
        if V1 in done_V1:
            print(f"=== Skipping V1={V1} (already in CSV) ===")
            continue

        print(f"\n=== Optimizing for V1={V1} ===")
        M_best = run_trials_with_refinement(t1=t1, t2=t2, V1=V1, V2=V2, **kwargs)
        G_best = compute_G(M_best)
        nn_bonds, nnn_bonds, nn_dirs, nnn_dirs = get_bond_lists(Lx, Ly)
        energy = energy_wick(G_best, nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2)
        staggered_order = compute_staggered_order(G_best, sublattices)

        key = f"t1_{t1:.3f}_V1_{V1:.3f}"

        # 安全写入 HDF5
        save_phi_to_h5_safe(h5_file, key, M_best)

        # 安全写入 CSV
        append_to_csv(csv_file, header, [t1, t2, V1, V2, energy.item(), staggered_order.item(), key])

        # 清理 JAX 缓存
        jax.clear_caches()

# -----------------------------------------
if __name__ == '__main__':
    start_time = time.time()

    sweep_v1(
        t1=1.0,
        t2=0.5,
        V2=0.0,
        V1_list=[0.28, 0.48, 0.60],
        Lx=4,
        Ly=4,
        rough_iter=300,
        refine_iter=1000,
        top_k=50,
        batch_size=16,
        lr=1e-2
    )

    end_time = time.time()
    elapsed = end_time - start_time
    print(f"\nTotal elapsed time: {elapsed/60:.2f} minutes ({elapsed:.1f} seconds)")
