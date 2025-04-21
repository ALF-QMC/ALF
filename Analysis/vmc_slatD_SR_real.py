#!/usr/bin/env python3
import jax
import jax.numpy as jnp
import numpy as np
from jax import vmap
from functools import partial
jax.config.update("jax_enable_x64", True)

def save_matrix_to_txt(matrix, filename="S_matrix.txt"):
    np_matrix = np.array(matrix)  # 将 JAX array 转为 NumPy array
    with open(filename, "w") as f:
        for row in np_matrix:
            row_str = "  ".join([f"{val:.3e}" for val in row])
            f.write(row_str + "\n")
    print(f"Matrix saved to {filename}")

def save_eigenvalues_to_txt(S, filename="S_eigvals.txt", n_per_line=1):
    eigvals = jnp.linalg.eigvalsh(S)
    eigvals_np = jnp.array(eigvals)  # 转换为 numpy-like array 便于保存
    with open(filename, "w") as f:
        for i in range(0, len(eigvals_np), n_per_line):
            line = "  ".join(f"{eigvals_np[j]:.10e}" for j in range(i, min(i + n_per_line, len(eigvals_np))))
            f.write(line + "\n")
    print(f"Eigenvalues written to {filename}")

# -----------------------------
# Configuration representation
# -----------------------------

def init_config(N, Nf, key):
    return jax.random.choice(key, N, (Nf,), replace=False)

# -----------------------------
# Utility for real/imag parameter handling
# -----------------------------

def reconstruct_M(params):
    return params["real"] + 1j * params["imag"]

# -----------------------------
# Local energy computation
# -----------------------------

def local_energy_from_W(xs, Ws, t_edges, V_edges):
    def energy_one(x, W):
        E_kin = 0.0
        for (i, j, t_ij) in t_edges:
            cond1 = (x == i).any() & ~(x == j).any()
            l1 = jnp.argmax(x == i)
            #kin1 = jnp.real(-t_ij * W[i, l1])
            kin1 = -t_ij * W[j, l1]
            ##jax.debug.print("Hop: x = {}, ", x.shape )
            ##jax.debug.print("Hop: i = {} to l1 = {} , kin = {}, ", i, l1, \
            ##        -t_ij * W[i, l1] )
            ##jax.debug.print("Hop: i = {} to l1 = {} , kin = {}, ", i, l1, \
            ##        -t_ij )
            E_kin += jnp.where(cond1, kin1, 0.0)
        
            cond2 = (x == j).any() & ~(x == i).any()
            l2 = jnp.argmax(x == j)
            #kin2 = jnp.real(-jnp.conj(t_ij) * W[j, l2])
            kin2 = -jnp.conj(t_ij) * W[i, l2]
            E_kin += jnp.where(cond2, kin2, 0.0)

        n = jnp.zeros(W.shape[0], dtype=jnp.float64)
        n = n.at[x].set(1.0)
        E_pot = 0.0
        for (i, j, V_ij) in V_edges:
            E_pot += V_ij * ((n[i] - 0.5) * (n[j] - 0.5))
            cond3 = (x == i).any() & (x == j).any()
            l1 = jnp.argmax(x == i)
            l2 = jnp.argmax(x == j)
            pot2 = -V_ij * W[j, l1]*W[i, l2]
            ##jax.debug.print("pot2: {}, W[i,i]: {}, W[j,j]: {}, W[i,j]: {}, W[j,i]: {}, cond3: {}", \
            ##        jnp.real(pot2), jnp.real(W[i,l1]), jnp.real(W[j,l2]), jnp.real(W[j,l1]), jnp.real(W[i,l2]), cond3 )
            E_pot += jnp.where(cond3, pot2, 0.0)
        
        #jax.debug.print("E: {}, ", E_kin + E_pot )

        return E_kin + E_pot

    ##jax.debug.print("Hop: xs = {}, ", xs.shape )
    energies = vmap(energy_one)(xs, Ws)
    return energies

# -----------------------------
# SR optimization step (real parameters)
# -----------------------------

def sr_update_real(params, xs, Ws, t_edges, V_edges, lr=0.05, damping=1e-4):
    M = reconstruct_M(params)
    D = M.shape[0] * M.shape[1]
    #jax.debug.print("grad_E = {}", grad_logs_combined)
    #jax.debug.print("M has NaN: {}", jnp.any(jnp.isnan(M)))
    #jax.debug.print("Number of NaNs in M: {}", jnp.sum(jnp.isnan(M)))
    #exit()

    # Flatten xs and Ws across time and walkers
    # Sampling stride
    stride = 50
    xs = xs[::stride]
    Ws = Ws[::stride]

    xs = xs.reshape(-1, xs.shape[-1])                      # (n_steps * n_walkers, Nf)
    Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])        # (n_steps * n_walkers, N, Nf)
    n_samples = xs.shape[0]

    def log_psi_grad(M, x):
        subM = M[x, :]
        subM_inv = jnp.linalg.inv(subM).T
        grad_log = jnp.zeros_like(M, dtype=M.dtype).at[x, :].set(subM_inv)
        #jax.debug.print("grad_log has NaN: {}", jnp.any(jnp.isnan(grad_log)))
        #jax.debug.print("Number of NaNs in grad_log: {}", jnp.sum(jnp.isnan(grad_log)))
        return grad_log

    ##def log_psi_grad(M, x):
    ##    subM = M[x, :]
    ##    subM_inv = jnp.linalg.inv(subM).T

    ##    def body(i, g):
    ##        return g.at[x[i], :].set(subM_inv[i, :])

    ##    grad_log = jax.lax.fori_loop(0, x.shape[0], body, jnp.zeros_like(M))
    ##    return grad_log
    
    grad_logs = vmap(log_psi_grad, in_axes=(None, 0))(M, xs)
    #jax.debug.print("grad_logs_combined has NaN: {}", jnp.any(jnp.isnan(grad_logs)))
    #jax.debug.print("Number of NaNs in grad_logs_combined: {}", jnp.sum(jnp.isnan(grad_logs)))
    #exit()
    grad_logs_flat = grad_logs.reshape((n_samples, -1))  # complex

    # Use correct directional derivatives:
    grad_logs_combined = jnp.concatenate([
        grad_logs_flat,
        1j * grad_logs_flat,
    ], axis=1)  # shape: (n_samples, 2D)
    
    E_loc = local_energy_from_W(xs, Ws, t_edges, V_edges)
    O_mean = jnp.mean(grad_logs_combined, axis=0, keepdims=True)
    #jax.debug.print("O_mean: {}", jnp.max(jnp.abs(O_mean)))
    E_mean = jnp.mean(E_loc)
    #jax.debug.print("E_mean: {}", E_mean)
    ElO = jnp.mean(jnp.conj(E_loc)[:, None] * grad_logs_combined, axis=0)
    #jax.debug.print("ElO: {}", jnp.max(jnp.abs(ElO)))

    grad_E = -2.0 * jnp.real(ElO - jnp.conj(E_mean) * O_mean.squeeze())
    #jax.debug.print("‖grad_E‖² = {}", jnp.real(jnp.vdot(grad_E, grad_E)))
    #jax.debug.print("grad_E: {}", jnp.max(jnp.abs(grad_E)))

    #jax.debug.print("grad_logs_combined.shape = {}", grad_logs_combined.shape)
    O_centered = grad_logs_combined - O_mean
    #jax.debug.print("O_mean.shape = {}", O_mean.shape)
    #jax.debug.print("O_centered.shape = {}", O_centered.shape)
    #jax.debug.print("O_centered max: {}", jnp.max(jnp.abs(O_centered)))

    #S = jnp.real((O_centered.T @ O_centered) / n_samples)
    S = jnp.real((jnp.conj(O_centered).T @ O_centered) / n_samples)
    #jax.debug.print("S.shape = {}", S.shape)

    # regularization
    # Preconditioned S
    diag_S = jnp.diag(S)
    scales = jnp.sqrt(diag_S)

    S_pc = S / (scales[:, None] * scales[None, :])
    S_pc += damping * jnp.eye(S.shape[0])

    # Preconditioned grad_E
    grad_E_pc = grad_E / scales

    #jax.debug.print("S shape: {}", S.shape)
    #jax.debug.print("S has NaN: {}", jnp.any(jnp.isnan(S)))
    #jax.debug.print("S max: {}, min: {}", jnp.max(S), jnp.min(S))
    #jax.debug.print("S trace: {}", jnp.trace(S))
    #save_matrix_to_txt(S[:D,:D], "S_matrix.txt")

    #ridge = damping * jnp.trace(S) / S.shape[0]
    #ridge = damping
    #S += ridge * jnp.eye(S.shape[0])

    #eigvals = jnp.linalg.eigvalsh(S[:D,:D])  # 只计算本征值（适合实对称矩阵）
    #jax.debug.print("S eigenvalues (min ~ max): {} ~ {}", jnp.min(eigvals), jnp.max(eigvals))
    #jax.debug.print("Smallest 5 eigvals: {}", eigvals[:5])
    #jax.debug.print("Largest 5 eigvals: {}", eigvals[-5:])
    #jax.debug.print("Smallest 100 eigvals: {}", eigvals[:100])
    #jax.debug.print("Largest 100 eigvals: {}", eigvals[-100:])
    #save_eigenvalues_to_txt(S[:D,:D])
    save_eigenvalues_to_txt(S_pc[:D,:D])

    #delta = jnp.linalg.solve(S, grad_E)
    delta_pc = jnp.linalg.solve(S_pc, grad_E_pc)
    # Map back to original space
    delta = delta_pc / scales
    #delta = jnp.linalg.solve(S[:D,:D], grad_E[:D])

    #D = M.shape[0] * M.shape[1]
    delta_real = delta[:D].reshape(M.shape)
    delta_imag = delta[D:].reshape(M.shape)
    #delta_real = grad_E[:D].reshape(M.shape)
    #delta_imag = grad_E[D:].reshape(M.shape)
    #jax.debug.print("delta_imag sum: {}", jnp.sum(delta_imag))

    #jax.debug.print("delta_real = {}", delta_real)

    #new_params = {
    #    "real": params["real"] - lr * delta_real,
    #    "imag": params["imag"] - lr * delta_imag,
    #}
    new_params = {
        "real": params["real"] + lr * delta_real,
        "imag": params["imag"] + lr * delta_imag,
    }
    return new_params

# -----------------------------
# Initial walker state
# -----------------------------

def init_walker_states(M, N, Nf, n_walkers, key):
    keys = jax.random.split(key, n_walkers)

    def update_best(state, trial):
        best_x, best_logabs, best_sign, best_W = state
        x_trial, logabs_trial, sign_trial, W_trial = trial

        cond = logabs_trial > best_logabs
        new_state = (
            x_trial,
            logabs_trial,
            sign_trial,
            W_trial,
        )
        return jax.lax.cond(cond, lambda _: new_state, lambda _: state, operand=None)

    def single_init(key):
        def body(i, state):
            key_i = jax.random.fold_in(key, i)
            x = init_config(N, Nf, key_i)
            subM = M[x, :]
            sign, logabs = jnp.linalg.slogdet(subM)
            W = jnp.linalg.solve(subM.T, M.T).T
            trial = (x, logabs, sign, W)
            return update_best(state, trial)

        # 初始值：最小 logabs
        x0 = init_config(N, Nf, key)
        subM0 = M[x0, :]
        sign0, logabs0 = jnp.linalg.slogdet(subM0)
        W0 = jnp.linalg.solve(subM0.T, M.T).T
        init_state = (x0, logabs0, sign0, W0)

        best_state = jax.lax.fori_loop(1, 100, body, init_state)
        return best_state

    ##def single_init(k):
    ##    x = init_config(N, Nf, k)
    ##    subM = M[x, :]
    ##    subM_inv = jnp.linalg.inv(subM)
    ##    sign, logabs = jnp.linalg.slogdet(subM)
    ##    W = jnp.linalg.solve(subM.T, M.T).T
    ##    
    ##    #jax.debug.print("sign = {}", sign)
    ##    #jax.debug.print("logabs = {}", logabs)
    ##    
    ##    #return x, logabs, sign, W
    ##    return x, logabs, sign, subM_inv

    walker_state = vmap(single_init)(keys)
    # 保存第一个 walker 的 W 矩阵到文本文件
    #W0 = walker_state[3][0]  # shape (N, Nf)
    W0 = walker_state[3][0]  # shape (N, Nf)
    W0_np = jnp.array(W0)    # 强制变成 numpy-like array 以便保存
    
    with open("W0_matrix.txt", "w") as f:
        for i in range(W0_np.shape[0]):
            row_str = "  ".join([f"{W0_np[i, j].real:.2f}" for j in range(W0_np.shape[1])])
            f.write(row_str + "\n")

    #jax.debug.print("x init : {}, ", walker_state[0])
    
    #print("W0 matrix saved to W0_matrix.txt")
    #exit()

    return vmap(single_init)(keys)

# -----------------------------
# Metropolis sampling with rank-1 update
# -----------------------------

def sample_batch_rank1_stable(M, walker_state, n_steps, key, reset_interval=50):
    Nf = M.shape[1]
    n_walkers = walker_state[0].shape[0]
    #jax.debug.print("M_shape: {}, ", M.shape)
    #exit()

    def single_step(carry, key_step):
        state, step = carry
        xs, logabs_vals, sign_vals, Ws = state
        #jax.debug.print("shape: xs={}, logabs_vals={}, sign_vals={}, Ws={}", xs.shape, \
        #        logabs_vals.shape, sign_vals.shape, Ws.shape)
        subkeys = jax.random.split(key_step, n_walkers)

        def single_update(x, logabs, sign, W, k):
            #jax.debug.print("shape: x={}, logabs_vals={}, sign_vals={}, W={}", x.shape, \
            #        logabs.shape, sign.shape, W.shape)

            key1, key2 = jax.random.split(k)
            idx = jax.random.randint(key1, (), 0, Nf)
            l1 = x[idx]

            occ_mask = jnp.zeros(M.shape[0], dtype=bool).at[x].set(True)
            unocc_mask = ~occ_mask
            weights = unocc_mask.astype(jnp.float64)
            probs = weights / jnp.sum(weights)
            k_site = jax.random.choice(key2, jnp.arange(M.shape[0]), p=probs)

            ratio = W[k_site, idx]
            abs_val = jnp.abs(ratio) + 1e-40
            log_ratio = jnp.log(abs_val)
            phase_ratio = ratio / abs_val

            delta = 2.0 * log_ratio
            accept = jax.random.uniform(key1) < jnp.exp(delta)
            #accept = False
            #jax.debug.print("delta : {}", jax.nn.one_hot(idx, Nf) )
            #jax.debug.print("W[k,:] : {}", W[k_site, :] )
            #jax.debug.print("W[k,:]-delta : {}", (W[k_site, :]-jax.nn.one_hot(idx, Nf)) )

            def accepted_update():
                x_new = x.at[idx].set(k_site)
                W_new = W - jnp.outer(W[:, idx], (W[k_site, :] - jax.nn.one_hot(idx, Nf)) / (ratio + 1e-40))
                #jax.debug.print("jnp.outer dtype: {}", (jnp.outer(W[:, idx], (W[k_site, :] - jax.nn.one_hot(l1, Nf)))).dtype )
                #jax.debug.print("NaN: W_new : {}, W : {}, ratio: {}, accept: {}", jnp.any(jnp.isnan(W_new)), \
                #        jnp.any(jnp.isnan(W)), ratio, accept )
                #jax.debug.print("idx: {}, x_idx: {}, k_site: {}", idx, l1, k_site )
                #jax.debug.print("max|W| = {}, max|W_new| = {}, D = {}, accept = {}", \
                #        jnp.max(jnp.abs(W)), jnp.max(jnp.abs(W_new)), \
                #        jnp.max(jnp.abs( W[k_site, :] - jax.nn.one_hot(idx, Nf) )), accept )
                W_new = W_new.astype(jnp.complex64)
                #jax.debug.print("W_new has NaN: {}, accept {}, ratio {}, ", jnp.any(jnp.isnan(W_new)), \
                #        accept, jnp.abs(ratio))
                return x_new, logabs + log_ratio, sign * phase_ratio, W_new

            def rejected_update():
                return x, logabs, sign, W

            return jax.lax.cond(accept, accepted_update, rejected_update)

        xs_new, logabs_new, sign_new, Ws_new = vmap(single_update)(xs, logabs_vals, sign_vals, Ws, subkeys)
        #jax.debug.print("xs_new : {}, ", xs_new)

        # Reset W every `reset_interval` steps
        def do_reset():
            def reset_W(x, W_in):
                #jax.debug.print("reset x shape: x={}", x.shape)
                #jax.debug.print("reset M shape: M={}", M.shape)
                subM = M[x, :]
                #subM_inv = jnp.linalg.inv(subM)
                sign_reset, logabs_reset = jnp.linalg.slogdet(subM)
                #jax.debug.print("logabs: reset = {}; in = {}", logabs_reset, logabs_new[0])
                #jax.debug.print("sign: reset = {}; in = {}", sign_reset, sign_new[0])
                #jax.debug.print("max: subM={}", jnp.max(subM))
                W_reset = jnp.linalg.solve(subM.T, M.T).T
                #W_reset = M@subM_inv
                #jax.debug.print("reset shape: W_reset={}", W_reset.shape)
                #jax.debug.print("W_reset has NaN : {}, ", jnp.any(jnp.isnan(W_reset)) )
                #W_test = W_reset[:,:] - W_in[:,:]
                #jax.debug.print("Diff W : {}, Max W_in: {}", jnp.sum(jnp.abs(W_test)), \
                #        jnp.max(jnp.abs(W_in)))
                return logabs_reset, sign_reset, W_reset
            logabs_reset, sign_reset, Ws_reset = vmap(reset_W)(xs_new, Ws_new)
            return xs_new, logabs_reset, sign_reset, Ws_reset

        def no_reset():
            return xs_new, logabs_new, sign_new, Ws_new

        #jax.debug.print("log: {}; sign = {}", logabs_new[0], sign_new[0])
        state_new = jax.lax.cond(step % reset_interval == 0, do_reset, no_reset)
        #jax.debug.print("log_reset: {}; sign_reset = {}", logabs_new[0], sign_new[0])

        return (state_new, step + 1), state_new

    keys = jax.random.split(key, n_steps)
    init_step = 0
    final_carry, hist = jax.lax.scan(single_step, (walker_state, init_step), keys)
    final_state = final_carry[0]

    #exit()
    return final_state, hist

#def sample_batch_rank1_stable(M, walker_state, n_steps, key):
#    Nf = M.shape[1]
#    n_walkers = walker_state[0].shape[0]
#
#    def single_step(state, key):
#        xs, logabs_vals, sign_vals, Ws = state
#        subkeys = jax.random.split(key, n_walkers)
#
#        def single_update(x, logabs, sign, W, k):
#            key1, key2 = jax.random.split(k)
#            idx = jax.random.randint(key1, (), 0, Nf)
#
#            occ_mask = jnp.zeros(M.shape[0], dtype=bool).at[x].set(True)
#            unocc_mask = ~occ_mask
#            weights = unocc_mask.astype(jnp.float64)
#            probs = weights / jnp.sum(weights)
#            k_site = jax.random.choice(key2, jnp.arange(M.shape[0]), p=probs)
#
#            ratio = W[k_site, idx]
#            #jax.debug.print("k_site = {}, idx = {}, ratio = {}", k_site, idx, ratio)
#            abs_val = jnp.abs(ratio) + 1e-40
#            log_ratio = jnp.log(abs_val)
#            phase_ratio = ratio / abs_val
#
#            delta = 2.0 * log_ratio
#            #jax.debug.print("exp(delta) = {}", jnp.exp(delta))
#            accept = jax.random.uniform(key1) < jnp.exp(delta)
#            #jax.debug.print("accept = {}, ratio = {}", accept, jnp.abs(ratio))
#
#            def accepted_update():
#                x_new = x.at[idx].set(k_site)
#                W_new = W - jnp.outer(W[:, idx], (W[k_site, :] - jax.nn.one_hot(idx, Nf)) / (ratio + 1e-40))
#                #jax.debug.print("ratio = {}", ratio)
#                W_new = W_new.astype(jnp.complex64)
#                jax.debug.print("W_new has NaN: {}, accept {}, ratio {}, ", jnp.any(jnp.isnan(W_new)), \
#                        accept, jnp.abs(ratio))
#                return x_new, logabs + log_ratio, sign * phase_ratio, W_new
#
#            def rejected_update():
#                return x, logabs, sign, W
#
#            return jax.lax.cond(accept, accepted_update, rejected_update)
#
#        new_state = vmap(single_update)(xs, logabs_vals, sign_vals, Ws, subkeys)
#        return new_state, new_state
#
#    keys = jax.random.split(key, n_steps)
#    final_state, hist = jax.lax.scan(single_step, walker_state, keys)
#    return final_state, hist

# -----------------------------
# Free Slater determinant initialization
# -----------------------------

def build_free_slater(Lx, Ly, t1, t2, m):
    N = 2 * Lx * Ly
    coords = []
    sublattices = []
    site_id = lambda x, y, s: 2 * (y % Ly * Lx + x % Lx) + s

    for y in range(Ly):
        for x in range(Lx):
            coords.append((x, y))
            coords.append((x + 0.5, y + 0.5))
            sublattices.extend([0, 1])

    H = jnp.zeros((N, N), dtype=jnp.complex64)
    for y in range(Ly):
        for x in range(Lx):
            ia = site_id(x, y, 0)
            ib = site_id(x, y, 1)
            H = H.at[ia, ib].set(-t1)
            H = H.at[ib, ia].set(-t1)
            H = H.at[ib, site_id(x+1, y, 0)].set(-t1)
            H = H.at[site_id(x+1, y, 0), ib].set(-t1)
            H = H.at[ib, site_id(x, y+1, 0)].set(-t1)
            H = H.at[site_id(x, y+1, 0), ib].set(-t1)
            H = H.at[ib, site_id(x+1, y+1, 0)].set(-t1)
            H = H.at[site_id(x+1, y+1, 0), ib].set(-t1)

            for s, sub in enumerate([0, 1]):
                i = site_id(x, y, sub)
                if sub == 0:
                    H = H.at[i, site_id(x+1, y, sub)].set(-t2)
                    H = H.at[site_id(x+1, y, sub), i].set(-t2)
                    H = H.at[i, site_id(x, y+1, sub)].set(t2)
                    H = H.at[site_id(x, y+1, sub), i].set(t2)
                else:
                    H = H.at[i, site_id(x+1, y, sub)].set(t2)
                    H = H.at[site_id(x+1, y, sub), i].set(t2)
                    H = H.at[i, site_id(x, y+1, sub)].set(-t2)
                    H = H.at[site_id(x, y+1, sub), i].set(-t2)

    sublattices = jnp.array(sublattices)
    for i in range(N):
        H = H.at[i, i].add(-m if sublattices[i] == 0 else m)

    H_np = jnp.array(H)
    with open("H_matrix.txt", "w") as f:
        for i in range(H_np.shape[0]):
            row_str = "  ".join([f"{H_np[i, j].real:.2f}" for j in range(H_np.shape[1])])
            f.write(row_str + "\n")
    print("Hamiltonian matrix written in matrix format to H_matrix.txt")

    evals, evecs = jnp.linalg.eigh(H)
    Nf = N // 2
    M0 = evecs[:, :Nf].astype(jnp.complex64)
    return M0

# -----------------------------
# Training loop
# -----------------------------

##def train(M_init, t_edges, V_edges, n_iter=100, n_walkers=64, n_steps=20, lr=0.05, n_thermal=50, key=jax.random.PRNGKey(42)):
##    N, Nf = M_init.shape
##    params = {
##        "real": M_init.real,
##        "imag": M_init.imag,
##    }
##
##    key, therm_key = jax.random.split(key)
##    #M = reconstruct_M(params)
##    #walker_state = init_walker_states(M, N, Nf, n_walkers, therm_key)
##
##    ##for _ in range(n_thermal):
##    ##    key, step_key = jax.random.split(key)
##    ##    M = reconstruct_M(params)
##    ##    walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, 100, step_key, reset_interval=50)
##    ##    #walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, 100, key=key, reset_interval=20)
##    ##    flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
##    ##    flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # shape (n_steps * n_walkers, N, Nf)
##    ##    energies = local_energy_from_W(flat_xs, flat_Ws, t_edges, V_edges)
##    ##    E_mean = jnp.mean(energies)
##    ##    
##    ##    print(f"Thermal: E = {E_mean:.6f}")
##    ##    exit()
##
##    ##for step in range(n_iter):
##    ##    key, step_key = jax.random.split(key)
##    ##    M = reconstruct_M(params)
##    ##    #walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, step_key)
##    ##    walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, key=key, reset_interval=50)
##
##    ##    params = sr_update_real(params, xs, Ws, t_edges, V_edges, lr=lr)
##
##    ##    flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
##    ##    flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # shape (n_steps * n_walkers, N, Nf)
##    ##    energies = local_energy_from_W(flat_xs, flat_Ws, t_edges, V_edges)
##    ##    E_mean = jnp.mean(energies)
##
##    ##    N = M.shape[0]  # 总格点数
##    ##    flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
##    ##    counts = jnp.sum(vmap(lambda x: jnp.bincount(x, length=N))(flat_xs), axis=0)
##    ##    N_tot = jnp.sum(counts) / (n_walkers * n_steps)
##
##    ##    print(f"Step {step:3d}: E = {E_mean:.6f}, N = {N_tot:.4f}")
##    
##    for step in range(n_iter):
##        M = reconstruct_M(params)
##        M = M.astype(jnp.complex64)
##        
##        #jax.debug.print("M imag : {}, ", M[0,:].imag)
##        jax.debug.print("M real : {}, ", M[0,:].real)
##
##        walker_state = init_walker_states(M, N, Nf, n_walkers, therm_key)
##    
##        ## warmup
##        walker_state, _ = sample_batch_rank1_stable(M, walker_state, n_thermal, key=therm_key, reset_interval=50)
##        
##        ## normal run
##        walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, key=key, reset_interval=50)
##
##        params = sr_update_real(params, xs, Ws, t_edges, V_edges, lr=lr)
##
##        flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
##        flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # shape (n_steps * n_walkers, N, Nf)
##        energies = local_energy_from_W(flat_xs, flat_Ws, t_edges, V_edges)
##        E_mean = jnp.mean(energies)
##
##        N = M.shape[0]  # 总格点数
##        flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
##        counts = jnp.sum(vmap(lambda x: jnp.bincount(x, length=N))(flat_xs), axis=0)
##        N_tot = jnp.sum(counts) / (n_walkers * n_steps)
##
##        print(f"Step {step:3d}: E = {E_mean:.6f}, N = {N_tot:.4f}")
##
##    return reconstruct_M(params)

###def train(M_init, t_edges, V_edges, n_iter=100, n_walkers=64, n_steps=20, lr=0.05, n_thermal=50, key=jax.random.PRNGKey(42)):
###    N, Nf = M_init.shape
###    params = {
###        "real": M_init.real,
###        "imag": M_init.imag,
###    }
###
###    # 构造 key tree
###    key, key_init, key_thermal, key_train = jax.random.split(key, 4)
###
###    # 初始化 walker
###    M = reconstruct_M(params).astype(jnp.complex64)
###    walker_state = init_walker_states(M, N, Nf, n_walkers, key_init)
###
###    # thermalization 步 key 拆分
###    thermal_keys = jax.random.split(key_thermal, n_thermal)
###    for i in range(n_thermal):
###        walker_state, _ = sample_batch_rank1_stable(M, walker_state, 1, key=thermal_keys[i], reset_interval=50)
###
###    # training 步 key 拆分
###    train_keys = jax.random.split(key_train, n_iter)
###    for step in range(n_iter):
###        M = reconstruct_M(params).astype(jnp.complex64)
###
###        # 每步 sampling 使用新的 key
###        walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, key=train_keys[step], reset_interval=50)
###
###        params = sr_update_real(params, xs, Ws, t_edges, V_edges, lr=lr)
###
###        flat_xs = xs.reshape(-1, xs.shape[-1])  # shape (n_steps * n_walkers, Nf)
###        flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # shape (n_steps * n_walkers, N, Nf)
###        energies = local_energy_from_W(flat_xs, flat_Ws, t_edges, V_edges)
###        E_mean = jnp.mean(energies)
###
###        counts = jnp.sum(vmap(lambda x: jnp.bincount(x, length=N))(flat_xs), axis=0)
###        N_tot = jnp.sum(counts) / (n_walkers * n_steps)
###
###        print(f"Step {step:3d}: E = {E_mean:.6f}, N = {N_tot:.4f}")
###
###    return reconstruct_M(params)

def train(M_init, t_edges, V_edges, n_iter=100, n_walkers=64, n_steps=20, lr=0.05, n_thermal=50, key=jax.random.PRNGKey(42)):
    N, Nf = M_init.shape
    params = {
        "real": M_init.real,
        "imag": M_init.imag,
    }

    # 预分配 train loop 的 key
    key, key_train = jax.random.split(key)
    train_keys = jax.random.split(key_train, n_iter)

    for step in range(n_iter):
        # 每一步分出 key_init, key_thermal, key_sample
        key_step = train_keys[step]
        key_init, key_thermal, key_sample = jax.random.split(key_step, 3)

        # 构造当前 Slater matrix
        M = reconstruct_M(params).astype(jnp.complex64)
        M, _ = jnp.linalg.qr(M)
        params = {
            "real": M.real,
            "imag": M.imag,
        }

        # 初始化 walker
        walker_state = init_walker_states(M, N, Nf, n_walkers, key_init)

        # thermalization
        walker_state, _ = sample_batch_rank1_stable(M, walker_state, n_thermal, key=key_thermal, reset_interval=50)

        # Metropolis sampling
        walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, key=key_sample, reset_interval=50)

        # SR 优化
        params = sr_update_real(params, xs, Ws, t_edges, V_edges, lr=lr)

        # 计算能量
        flat_xs = xs.reshape(-1, xs.shape[-1])      # (n_steps * n_walkers, Nf)
        flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # (n_steps * n_walkers, N, Nf)
        energies = local_energy_from_W(flat_xs, flat_Ws, t_edges, V_edges)
        E_mean = jnp.mean(energies)

        # 粒子数统计
        counts = jnp.sum(vmap(lambda x: jnp.bincount(x, length=N))(flat_xs), axis=0)
        N_tot = jnp.sum(counts) / (n_walkers * n_steps)

        print(f"Step {step:3d}: E = {E_mean:.6f}, N = {N_tot:.4f}")

    return reconstruct_M(params)

# -----------------------------
# Main function
# -----------------------------

def main():
    Lx, Ly = 4, 4
    t1, t2, V1, V2 = 1.0, 0.5, 2.0, 0.0
    m = 0.01
    key = jax.random.PRNGKey(0)
    N = 2 * Lx * Ly

    def site_id(x, y, s):
        return 2 * ((y % Ly) * Lx + (x % Lx)) + s

    t_edges, V_edges = [], []
    for y in range(Ly):
        for x in range(Lx):
            ia = site_id(x, y, 0)
            ib = site_id(x, y, 1)
            t_edges.extend([
                (ib, ia, t1),
                (ib, site_id(x+1, y, 0), t1),
                (ib, site_id(x, y+1, 0), t1),
                (ib, site_id(x+1, y+1, 0), t1),
                (site_id(x, y, 0), site_id(x+1, y, 0),  t2),
                (site_id(x, y, 0), site_id(x, y+1, 0), -t2),
                (site_id(x, y, 1), site_id(x+1, y, 1), -t2),
                (site_id(x, y, 1), site_id(x, y+1, 1),  t2),
            ])
            V_edges.extend([
                (ib, ia, V1),
                (ib, site_id(x+1, y, 0), V1),
                (ib, site_id(x, y+1, 0), V1),
                (ib, site_id(x+1, y+1, 0), V1),
                (site_id(x, y, 0), site_id(x+1, y, 0), V2),
                (site_id(x, y, 0), site_id(x, y+1, 0), V2),
                (site_id(x, y, 1), site_id(x+1, y, 1), V2),
                (site_id(x, y, 1), site_id(x, y+1, 1), V2),
            ])

    M0 = build_free_slater(Lx, Ly, t1, t2, m)
    
    M_final = train(M0, t_edges, V_edges, n_iter=20, n_walkers=1000, n_steps=2000, lr=1e-2, n_thermal=100, key=key)
    return M_final

if __name__ == "__main__":
    main()
