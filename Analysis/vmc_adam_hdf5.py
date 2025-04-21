#!/usr/bin/env python3
import jax
import jax.numpy as jnp
import numpy as np
import optax
import h5py
from jax.scipy.linalg import eigh
jax.config.update("jax_enable_x64", True)

def save_phi_trial(M, filename="phi_trial.h5"):
    """
    保存 Slater 行列式 M 为 (2, ndim, n_part) 的 float64 数组，确保虚部不为零。
    使用 h5dump 检查时也能看到。
    """
    M = np.asarray(M, dtype=np.complex128)  # 保证是 numpy complex128 类型
    ndim, n_part = M.shape

    # 构造 (2, ndim, n_part) 的 float64 数组，第一维 0: 实部, 1: 虚部
    data = np.stack([M.real, M.imag], axis=-1)

    # 写入 HDF5 文件
    with h5py.File(filename, "w") as f:
        f.create_dataset("phi_trial", data=data, dtype=np.float64)

# -----------------------------------------
# Checkerboard Lattice Construction
# -----------------------------------------
def make_checkerboard_lattice(Lx, Ly):
    N_sites = 2 * Lx * Ly  # two sublattices per unit cell
    coords = []
    sublattices = []
    for y in range(Ly):
        for x in range(Lx):
            # a sublattice at (x, y)
            coords.append((x, y))
            sublattices.append(0)
            # b sublattice at (x+0.5, y+0.5)
            coords.append((x + 0.5, y + 0.5))
            sublattices.append(1)
    coords = jnp.array(coords)
    sublattices = jnp.array(sublattices)
    return N_sites, coords, sublattices

# -----------------------------------------
# Hopping and Interaction Bond List
# -----------------------------------------
def get_bond_lists(Lx, Ly):
    def site_index(x, y, sub):
        return 2 * ((y % Ly) * Lx + (x % Lx)) + sub

    nn_bonds = []
    nn_dirs = []
    nnn_bonds = []
    nnn_dirs = []

    for y in range(Ly):
        for x in range(Lx):
            ia = site_index(x, y, 0)
            ib = site_index(x, y, 1)

            # nearest neighbors (a-b): 4 directions
            nn_bonds.append((ib, ia))
            nn_dirs.append(( 1, 'n1'))
            nn_bonds.append((ib, site_index(x, y+1, 0)))
            nn_dirs.append((-1, 'n1'))
            nn_bonds.append((ib, site_index(x+1, y, 0)))
            nn_dirs.append((-1, 'n1'))
            nn_bonds.append((ib, site_index(x+1, y+1, 0)))
            nn_dirs.append(( 1, 'n1'))

            # next-nearest neighbors (within same sublattice)
            # a sublattice
            ja_x = site_index(x+1, y, 0)
            ja_y = site_index(x, y+1, 0)
            nnn_bonds.append((ia, ja_x))
            nnn_dirs.append((0, 'x'))
            nnn_bonds.append((ia, ja_y))
            nnn_dirs.append((0, 'y'))
            # b sublattice
            jb_x = site_index(x+1, y, 1)
            jb_y = site_index(x, y+1, 1)
            nnn_bonds.append((ib, jb_x))
            nnn_dirs.append((1, 'x'))
            nnn_bonds.append((ib, jb_y))
            nnn_dirs.append((1, 'y'))

    return list(set(nn_bonds)), nnn_bonds, nn_dirs, nnn_dirs

# -----------------------------------------
# Free Hamiltonian Construction
# -----------------------------------------
def make_free_hamiltonian(N, nn_bonds, nnn_bonds, nn_dirs, nnn_dirs, sublattices, t1, t2, m):
    H = jnp.zeros((N, N), dtype=jnp.complex64)
    #for i, j in nn_bonds:
    #    H = H.at[i, j].set(-t1)
    #    H = H.at[j, i].set(-t1)
    for (i, j), (sub, direction) in zip(nn_bonds, nn_dirs):
        phase = -t1*jnp.exp(sub*jnp.pi/4.0)
        H = H.at[i, j].set(phase)
        H = H.at[j, i].set(jnp.conj(phase))
    for (i, j), (sub, direction) in zip(nnn_bonds, nnn_dirs):
        sign = -t2 if (sub == 0 and direction == 'x') or (sub == 1 and direction == 'y') else t2
        H = H.at[i, j].set(sign)
        H = H.at[j, i].set(sign)
    for i in range(N):
        H = H.at[i, i].add((-m if sublattices[i] == 0 else m))
    return H

# -----------------------------------------
# Slater Green Function and Wick Energy
# -----------------------------------------
def compute_G(M):
    Q, _ = jnp.linalg.qr(M)
    return Q @ Q.conj().T

def energy_wick(G, nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2):
    E_kin = 0.0
    for i, j in nn_bonds:
        E_kin += -t1 * (G[i, j] + G[j, i])
    for (i, j), (sub, direction) in zip(nnn_bonds, nnn_dirs):
        sign = -t2 if (sub == 0 and direction == 'x') or (sub == 1 and direction == 'y') else t2
        E_kin += sign * (G[i, j] + G[j, i])

    E_int = 0.0
    for (i, j), V in zip(nn_bonds + nnn_bonds, [V1]*len(nn_bonds) + [V2]*len(nnn_bonds)):
        ni = G[i, i]
        nj = G[j, j]
        E_int += V * ((ni - 0.5)*(nj - 0.5) - G[i,j]*G[j,i])

    return jnp.real(E_kin + E_int)

# -----------------------------------------
# Main Optimization Loop
# -----------------------------------------
def run_vmc(Lx=4, Ly=4, t1=1.0, t2=0.5, V1=1.40, V2=0.00, m=0.2, lr=1e-2, n_iter=600):
    N, coords, sublattices = make_checkerboard_lattice(Lx, Ly)
    Nf = N // 2
    #nn_bonds, nnn_bonds, nnn_dirs = get_bond_lists(Lx, Ly)
    nn_bonds, nnn_bonds, nn_dirs, nnn_dirs = get_bond_lists(Lx, Ly)
    
    #H0 = make_free_hamiltonian(N, nn_bonds, nnn_bonds, nnn_dirs, sublattices, t1, t2, m)
    H0 = make_free_hamiltonian(N, nn_bonds, nnn_bonds, nn_dirs, nnn_dirs, sublattices, t1, t2, m)
    
    evals, evecs = eigh(H0)
    M0 = evecs[:, :Nf]

    ## random init
    np.random.seed(2143)  # 设置种子保证可重复
    M_real = np.random.randn(N, Nf)
    M_imag = np.random.randn(N, Nf)
    M0 = M_real + 1j * M_imag

    def pack(M):
        return jnp.concatenate([M.real.flatten(), M.imag.flatten()])

    def unpack(p):
        M_real = p[:N*Nf].reshape((N, Nf))
        M_imag = p[N*Nf:].reshape((N, Nf))
        return M_real + 1j * M_imag

    #ap = np.loadtxt('params.txt')
    #params = ap
    #M0 = unpack(params)

    params = pack(M0)
    opt = optax.adam(lr)
    #opt = optax.adamw(learning_rate=1e-2, weight_decay=1e-4)
    opt_state = opt.init(params)

    @jax.jit
    def step(params, opt_state):
        M = unpack(params)
        G = compute_G(M)
        loss = energy_wick(G, nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2)
        grads = jax.grad(lambda p: energy_wick(compute_G(unpack(p)), nn_bonds, nnn_bonds, nnn_dirs, V1, V2, t1, t2))(params)
        updates, opt_state = opt.update(grads, opt_state) ## adam
        #updates, opt_state = opt.update(grads, opt_state, params=params) ## adamw
        params = optax.apply_updates(params, updates)
        return params, opt_state, loss

    for i in range(n_iter):
        params, opt_state, loss = step(params, opt_state)
        if i % 1 == 0:
            G_diag = jnp.real(jnp.diag(compute_G(unpack(params))))
            N_total = jnp.sum(G_diag)
            print(f"Step {i}, Energy = {loss:.6f}, N = {N_total:.4f}")

    M_opt = unpack(params)
    #np.savetxt('params.txt', params)
    
    import matplotlib.pyplot as plt

    # 画出粒子密度分布图
    G_opt = compute_G(M_opt)
    density = jnp.real(jnp.diag(G_opt))
    x = coords[:, 0]
    y = coords[:, 1]

    plt.figure(figsize=(6, 5))
    sc = plt.scatter(x, y, c=density, cmap="viridis", s=100, edgecolors='k')
    plt.colorbar(sc, label="Density $n_i$")
    plt.title("Real-space density distribution")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis('equal')
    plt.tight_layout()
    plt.show()
    
    save_phi_trial(M_opt, filename="phi_trial.h5")

    return M_opt

if __name__ == '__main__':
    M_final = run_vmc()
