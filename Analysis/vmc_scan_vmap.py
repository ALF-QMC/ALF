#!/usr/bin/env python3

import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import os
from functools import partial

jax.config.update("jax_enable_x64", True)

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
    return jnp.array(rows), jnp.array(cols), jnp.array(values)

def build_bond_mapping(Lx, Ly, pbc=True):
    bond_list = []
    for ix in range(Lx):
        for iy in range(Ly):
            jb = get_index(ix, iy, 1, Lx, Ly)
            offsets = [(0, 0), (0, 1), (1, 0), (1, 1)]
            for a, (dx, dy) in enumerate(offsets):
                ia_x = (ix + dx) % Lx if pbc else ix + dx
                ia_y = (iy + dy) % Ly if pbc else iy + dy
                if 0 <= ia_x < Lx and 0 <= ia_y < Ly:
                    ia = get_index(ia_x, ia_y, 0, Lx, Ly)
                    bond_list.append((ia, jb))
    return jnp.array(bond_list, dtype=jnp.int32)

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
                    ib1 = get_index(ix, iy, 1, Lx, Ly)
                    ib2 = get_index(ix2, iy2, 1, Lx, Ly)
                    bond_list.append((ib1, ib2))
    return jnp.array(bond_list, dtype=jnp.int32)

@partial(jax.jit, static_argnums=(9,))
def build_hamiltonian(phi1, phi2, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite):
    H = jnp.zeros((Nsite, Nsite))
    H = H.at[rows, cols].set(values)
    for k, (ia, ib) in enumerate(bond_pairs1):
        H = H.at[ia, ia].add(-V1 * phi1[k])
        H = H.at[ib, ib].add(+V1 * phi1[k])
    for k, (i1, i2) in enumerate(bond_pairs2):
        H = H.at[i1, i1].add(-V2 * phi2[k])
        H = H.at[i2, i2].add(+V2 * phi2[k])
    H = H + 1e-8 * jnp.eye(Nsite)
    return H

def total_density_matrix(H, Nelec):
    eigvals, eigvecs = jnp.linalg.eigh(H)
    F_occ = jnp.where(jnp.arange(H.shape[0]) < Nelec, 1.0, 0.0)
    n_matrix = eigvecs @ jnp.diag(F_occ) @ eigvecs.T.conj()
    return n_matrix

@jax.jit
def compute_phi_from_density(n_in, bond_pairs1, bond_pairs2):
    phi1 = n_in[bond_pairs1[:, 0]] - n_in[bond_pairs1[:, 1]]
    phi2 = n_in[bond_pairs2[:, 0]] - n_in[bond_pairs2[:, 1]]
    return phi1, phi2

@partial(jax.jit, static_argnums=(7,))
def loss_fn(V, n_in, rows, cols, values, bond_pairs1, bond_pairs2, Nsite, Nelec):
    V1, V2 = V[0], V[1]
    phi1_in, phi2_in = compute_phi_from_density(n_in, bond_pairs1, bond_pairs2)
    H = build_hamiltonian(phi1_in, phi2_in, rows, cols, values, bond_pairs1, bond_pairs2, V1, V2, Nsite)
    n_matrix = total_density_matrix(H, Nelec)
    n_out = jnp.real(jnp.diag(n_matrix))
    phi1_out = n_out[bond_pairs1[:, 0]] - n_out[bond_pairs1[:, 1]]
    phi2_out = n_out[bond_pairs2[:, 0]] - n_out[bond_pairs2[:, 1]]
    return jnp.mean((phi1_out - phi1_in) ** 2) + jnp.mean((phi2_out - phi2_in) ** 2)

def grid_search_V(n_in, Lx, Ly, t1, t2p, Nelec, V1_range, V2_range, dV1, dV2, pbc=True):
    Nsite = int(2 * Lx * Ly)
    bond_pairs1 = build_bond_mapping(Lx, Ly, pbc)
    bond_pairs2 = build_bond_mapping_v2(Lx, Ly, pbc)
    rows, cols, values = build_checkerboard_neighbors(Lx, Ly, t1, t2p, pbc)

    # 1) 构造 1D 网格
    V1_vals = jnp.arange(V1_range[0], V1_range[1] + dV1, dV1)
    V2_vals = jnp.arange(V2_range[0], V2_range[1] + dV2, dV2)

    # 2) 构造 2D 网格 (len(V1_vals), len(V2_vals)) 并 flatten
    V1_mesh, V2_mesh = jnp.meshgrid(V1_vals, V2_vals, indexing='ij')
    V_grid = jnp.stack([V1_mesh.ravel(), V2_mesh.ravel()], axis=1)  # shape (N, 2)

    # 3) 单点 loss 函数，jit 编译
    @jax.jit
    def single_loss(V):
        return loss_fn(V, n_in, rows, cols, values, bond_pairs1, bond_pairs2, Nsite, Nelec)

    # 4) 并行评估多个 (V1, V2)
    multi_loss_fn = jax.vmap(single_loss, in_axes=(0,))
    losses = multi_loss_fn(V_grid)  # shape (N,)

    # 5) 将结果 reshape 回 2D
    loss_matrix = losses.reshape((len(V1_vals), len(V2_vals)))

    # 6) 找到最优 index
    idx_flat = jnp.argmin(losses)
    best_loss = losses[idx_flat]
    best_i = idx_flat // len(V2_vals)
    best_j = idx_flat % len(V2_vals)
    best_V1 = V1_vals[best_i]
    best_V2 = V2_vals[best_j]

    print(f"\nBest V: V1 = {best_V1:.3f}, V2 = {best_V2:.3f}, Loss = {best_loss:.6e}")
    end = time.perf_counter()
    print(f"Total runtime = {end - start:.4f} s")

    # 7) 画图保存
    import os
    os.makedirs("results", exist_ok=True)

    plt.figure(figsize=(6, 5))
    plt.imshow(jnp.log(loss_matrix.T), origin='lower',
               extent=[V1_vals[0], V1_vals[-1], V2_vals[0], V2_vals[-1]],
               aspect='auto', cmap='viridis')
    plt.colorbar(label='Loss')
    plt.xlabel("V1")
    plt.ylabel("V2")
    plt.title("Loss Landscape [vmap]")
    plt.tight_layout()
    plt.savefig("results/loss_heatmap.png", dpi=150)
    plt.show()

    return float(best_V1), float(best_V2)


if __name__ == "__main__":
    import time
    start = time.perf_counter()
    filename = "density.dat"
    n_in = jnp.array(np.loadtxt(filename))

    best_V = grid_search_V(
        n_in=n_in,
        Lx=12, Ly=12,
        t1=1.0, t2p=0.5,
        Nelec=16,
        V1_range=(0.0, 10.0),
        V2_range=(0.0, 10.0),
        dV1=0.1, dV2=0.1,
        pbc=True
    )
