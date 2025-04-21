#!/usr/bin/env python3
import jax
import jax.numpy as jnp
import numpy as np
from jax import vmap
from functools import partial
jax.config.update("jax_enable_x64", True)

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
@jax.jit
def local_energy_from_W(xs, Ws, t_edges, V_edges):
    def energy_one(x, W):
        n = jnp.zeros(W.shape[0])
        n = n.at[x].set(1.0)

        def kin_term(edge):
            i, j, t = edge
            i = i.astype(jnp.int32)
            j = j.astype(jnp.int32)
            cond1 = (x == i).any() & ~(x == j).any()
            cond2 = (x == j).any() & ~(x == i).any()
            l1 = jnp.argmax(x == i)
            l2 = jnp.argmax(x == j)
            e1 = -t * W[j, l1]
            e2 = -jnp.conj(t) * W[i, l2]
            return jnp.where(cond1, e1, 0.0) + jnp.where(cond2, e2, 0.0)

        def pot_term(edge):
            i, j, V = edge
            i = i.astype(jnp.int32)
            j = j.astype(jnp.int32)
            e1 = V * ((n[i] - 0.5) * (n[j] - 0.5))
            cond = (x == i).any() & (x == j).any()
            l1 = jnp.argmax(x == i)
            l2 = jnp.argmax(x == j)
            e2 = -V * W[j, l1] * W[i, l2]
            return e1 + jnp.where(cond, e2, 0.0)

        E_kin = jnp.sum(jax.vmap(kin_term)(t_edges))
        E_pot = jnp.sum(jax.vmap(pot_term)(V_edges))
        return E_kin + E_pot

    return jax.vmap(energy_one)(xs, Ws)

# -----------------------------
# SR optimization step (real parameters)
# -----------------------------
@jax.jit
def batched_log_psi_grad(M, xs):
    def log_psi_grad(M, x):
        subM = M[x, :]
        subM_inv = jnp.linalg.inv(subM).T
        return jnp.zeros_like(M, dtype=M.dtype).at[x, :].set(subM_inv)
    return vmap(log_psi_grad, in_axes=(None, 0))(M, xs)

#@jax.jit
#def sr_update_real(params, xs, Ws, t_edges, V_edges, lr=0.05, damping=1e-4):
#    M = reconstruct_M(params)
#    D = M.shape[0] * M.shape[1]
#
#    # Sampling stride
#    stride = 10
#    xs = xs[::stride]
#    Ws = Ws[::stride]
#
#    xs = xs.reshape(-1, xs.shape[-1])
#    Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])
#    n_samples = xs.shape[0]
#
#    grad_logs = batched_log_psi_grad(M, xs)
#    grad_logs_flat = grad_logs.reshape((n_samples, -1))  # shape: (n_samples, D)
#
#    grad_logs_combined = jnp.concatenate([
#        grad_logs_flat,
#        1j * grad_logs_flat,
#    ], axis=1)  # shape: (n_samples, 2D)
#
#    E_loc = local_energy_from_W(xs, Ws, t_edges, V_edges)
#    O_mean = jnp.mean(grad_logs_combined, axis=0, keepdims=True)
#    E_mean = jnp.mean(E_loc)
#    ElO = jnp.mean(jnp.conj(E_loc)[:, None] * grad_logs_combined, axis=0)
#
#    grad_E = -2.0 * jnp.real(ElO - jnp.conj(E_mean) * O_mean.squeeze())
#
#    O_centered = grad_logs_combined - O_mean
#    S = jnp.real((jnp.conj(O_centered).T @ O_centered) / n_samples)
#
#    # Preconditioned S
#    #diag_S = jnp.clip(jnp.diag(S), a_min=1e-8)  # avoid divide-by-zero
#    diag_S = jnp.diag(S)  # avoid divide-by-zero
#    scales = jnp.sqrt(diag_S)
#    S_pc = S / (scales[:, None] * scales[None, :])
#    S_pc += damping * jnp.eye(S.shape[0])
#
#    grad_E_pc = grad_E / scales
#    delta_pc = jnp.linalg.solve(S_pc, grad_E_pc)
#    delta = delta_pc / scales
#
#    delta_real = delta[:D].reshape(M.shape)
#    delta_imag = delta[D:].reshape(M.shape)
#
#    return {
#        "real": params["real"] + lr * delta_real,
#        "imag": params["imag"] + lr * delta_imag,
#    }

### -----------------------------
### SR optimization step (real parameters)
### -----------------------------
##
##def sr_update_real(params, xs, Ws, t_edges, V_edges, lr=0.05, damping=1e-4, batch_size=1000):
##    M = reconstruct_M(params)
##    D = M.shape[0] * M.shape[1]
##
##    # Sampling stride
##    stride = 10
##    xs = xs[::stride]
##    Ws = Ws[::stride]
##
##    xs = xs.reshape(-1, xs.shape[-1])
##    Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])
##    n_samples = xs.shape[0]
##
##    S_accum = jnp.zeros((2 * D, 2 * D))
##    grad_E_accum = jnp.zeros((2 * D,))
##    count = 0
##
##    for i in range(0, n_samples, batch_size):
##        xs_batch = xs[i:i+batch_size]
##        Ws_batch = Ws[i:i+batch_size]
##        n_b = xs_batch.shape[0]
##
##        grad_logs = batched_log_psi_grad(M, xs_batch)
##        grad_logs_flat = grad_logs.reshape((n_b, -1))
##
##        grad_logs_combined = jnp.concatenate([
##            grad_logs_flat,
##            1j * grad_logs_flat,
##        ], axis=1)  # shape: (n_b, 2D)
##
##        E_loc = local_energy_from_W(xs_batch, Ws_batch, t_edges, V_edges)
##        O_mean = jnp.mean(grad_logs_combined, axis=0, keepdims=True)
##        E_mean = jnp.mean(E_loc)
##        ElO = jnp.mean(jnp.conj(E_loc)[:, None] * grad_logs_combined, axis=0)
##
##        grad_E = -2.0 * jnp.real(ElO - jnp.conj(E_mean) * O_mean.squeeze())
##        O_centered = grad_logs_combined - O_mean
##
##        S_b = jnp.real((jnp.conj(O_centered).T @ O_centered) / n_b)
##
##        S_accum += S_b * n_b
##        grad_E_accum += grad_E * n_b
##        count += n_b
##
##    S = S_accum / count
##    grad_E = grad_E_accum / count
##
##    diag_S = jnp.clip(jnp.diag(S), a_min=1e-10)
##    scales = jnp.sqrt(diag_S)
##    S_pc = S / (scales[:, None] * scales[None, :])
##    S_pc += damping * jnp.eye(S.shape[0])
##
##    grad_E_pc = grad_E / scales
##    delta_pc = jnp.linalg.solve(S_pc, grad_E_pc)
##    delta = delta_pc / scales
##
##    delta_real = delta[:D].reshape(M.shape)
##    delta_imag = delta[D:].reshape(M.shape)
##
##    return {
##        "real": params["real"] + lr * delta_real,
##        "imag": params["imag"] + lr * delta_imag,
##    }

##def sr_update_real_minibatch(params, xs, Ws, t_edges, V_edges, batch_size=1000, lr=0.05, damping=1e-4):
##    M = reconstruct_M(params)
##    D = M.shape[0] * M.shape[1]
##
##    xs = xs.reshape(-1, xs.shape[-1])
##    Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])
##    n_samples = xs.shape[0]
##    n_batches = n_samples // batch_size
##
##    def batch_fn(carry, idx):
##        grad_E_sum, S_sum = carry
##        start = idx * batch_size
##        xs_b = jax.lax.dynamic_slice(xs, (start, 0), (batch_size, xs.shape[1]))
##        Ws_b = jax.lax.dynamic_slice(Ws, (start, 0, 0), (batch_size, Ws.shape[1], Ws.shape[2]))
##
##        grad_logs = batched_log_psi_grad(M, xs_b)
##        grad_logs_flat = grad_logs.reshape((batch_size, -1))
##        grad_logs_combined = jnp.concatenate([grad_logs_flat, 1j * grad_logs_flat], axis=1)
##
##        E_loc = local_energy_from_W(xs_b, Ws_b, t_edges, V_edges)
##        O_mean = jnp.mean(grad_logs_combined, axis=0, keepdims=True)
##        E_mean = jnp.mean(E_loc)
##        ElO = jnp.mean(jnp.conj(E_loc)[:, None] * grad_logs_combined, axis=0)
##
##        grad_E = -2.0 * jnp.real(ElO - jnp.conj(E_mean) * O_mean.squeeze())
##        O_centered = grad_logs_combined - O_mean
##        S = jnp.real((jnp.conj(O_centered).T @ O_centered) / batch_size)
##
##        return (grad_E_sum + grad_E, S_sum + S), None
##
##    (grad_E_total, S_total), _ = jax.lax.scan(batch_fn,
##                                              (jnp.zeros(2 * D), jnp.zeros((2 * D, 2 * D))),
##                                              jnp.arange(n_batches))
##    grad_E_total /= n_batches
##    S_total /= n_batches
##
##    diag_S = jnp.clip(jnp.diag(S_total), a_min=1e-10)
##    scales = jnp.sqrt(diag_S)
##    S_pc = S_total / (scales[:, None] * scales[None, :])
##    S_pc += damping * jnp.eye(S_pc.shape[0])
##
##    grad_E_pc = grad_E_total / scales
##    delta_pc = jnp.linalg.solve(S_pc, grad_E_pc)
##    delta = delta_pc / scales
##
##    delta_real = delta[:D].reshape(M.shape)
##    delta_imag = delta[D:].reshape(M.shape)
##
##    return {
##        "real": params["real"] + lr * delta_real,
##        "imag": params["imag"] + lr * delta_imag,
##    }
##
##    from jax import jit

@partial(jax.jit, static_argnames=['batch_size'])
def sr_update_real_minibatch(params, xs, Ws, t_edges, V_edges, batch_size=1000, lr=0.05, damping=1e-4):
    M = reconstruct_M(params)
    D = M.shape[0] * M.shape[1]

        # Apply stride
    stride = 50
    xs = xs[::stride]
    Ws = Ws[::stride]

    xs = xs.reshape(-1, xs.shape[-1])
    Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])
    n_samples = xs.shape[0]
    n_batches = n_samples // batch_size

    def batch_fn(carry, idx):
        grad_E_sum, S_sum = carry
        start = idx * batch_size
        xs_b = jax.lax.dynamic_slice(xs, (start, 0), (batch_size, xs.shape[1]))
        Ws_b = jax.lax.dynamic_slice(Ws, (start, 0, 0), (batch_size, Ws.shape[1], Ws.shape[2]))

        grad_logs = batched_log_psi_grad(M, xs_b)
        grad_logs_flat = grad_logs.reshape((batch_size, -1))
        grad_logs_combined = jnp.concatenate([grad_logs_flat, 1j * grad_logs_flat], axis=1)

        E_loc = local_energy_from_W(xs_b, Ws_b, t_edges, V_edges)
        O_mean = jnp.mean(grad_logs_combined, axis=0, keepdims=True)
        E_mean = jnp.mean(E_loc)
        ElO = jnp.mean(jnp.conj(E_loc)[:, None] * grad_logs_combined, axis=0)

        grad_E = -2.0 * jnp.real(ElO - jnp.conj(E_mean) * O_mean.squeeze())
        O_centered = grad_logs_combined - O_mean
        S = jnp.real((jnp.conj(O_centered).T @ O_centered) / batch_size)

        return (grad_E_sum + grad_E, S_sum + S), None

    (grad_E_total, S_total), _ = jax.lax.scan(
        batch_fn,
        (jnp.zeros(2 * D, dtype=jnp.float64), jnp.zeros((2 * D, 2 * D), dtype=jnp.float64)),
        jnp.arange(n_batches)
    )

    grad_E_total /= n_batches
    S_total /= n_batches

    diag_S = jnp.clip(jnp.diag(S_total), a_min=1e-10)
    scales = jnp.sqrt(diag_S)
    S_pc = S_total / (scales[:, None] * scales[None, :])
    S_pc += damping * jnp.eye(S_pc.shape[0])

    grad_E_pc = grad_E_total / scales
    delta_pc = jnp.linalg.solve(S_pc, grad_E_pc)
    delta = delta_pc / scales

    delta_real = delta[:D].reshape(M.shape)
    delta_imag = delta[D:].reshape(M.shape)

    return {
        "real": params["real"] + lr * delta_real,
        "imag": params["imag"] + lr * delta_imag,
    }

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

    walker_state = vmap(single_init)(keys)
    # 保存第一个 walker 的 W 矩阵到文本文件
    W0 = walker_state[3][0]  # shape (N, Nf)
    W0_np = jnp.array(W0)    # 强制变成 numpy-like array 以便保存
    
    with open("W0_matrix.txt", "w") as f:
        for i in range(W0_np.shape[0]):
            row_str = "  ".join([f"{W0_np[i, j].real:.6f}" for j in range(W0_np.shape[1])])
            f.write(row_str + "\n")

    return vmap(single_init)(keys)

# -----------------------------
# Metropolis sampling with rank-1 update
# -----------------------------

def sample_batch_rank1_stable(M, walker_state, n_steps, key, reset_interval=50):
    Nf = M.shape[1]
    n_walkers = walker_state[0].shape[0]

    def single_step(carry, key_step):
        state, step = carry
        xs, logabs_vals, sign_vals, Ws = state
        subkeys = jax.random.split(key_step, n_walkers)

        def single_update(x, logabs, sign, W, k):

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

            def accepted_update():
                x_new = x.at[idx].set(k_site)
                W_new = W - jnp.outer(W[:, idx], (W[k_site, :] - jax.nn.one_hot(idx, Nf)) / (ratio + 1e-40))
                W_new = W_new.astype(jnp.complex64)
                return x_new, logabs + log_ratio, sign * phase_ratio, W_new

            def rejected_update():
                return x, logabs, sign, W

            return jax.lax.cond(accept, accepted_update, rejected_update)

        xs_new, logabs_new, sign_new, Ws_new = vmap(single_update)(xs, logabs_vals, sign_vals, Ws, subkeys)

        # Reset W every `reset_interval` steps
        def do_reset():
            def reset_W(x, W_in):
                subM = M[x, :]
                sign_reset, logabs_reset = jnp.linalg.slogdet(subM)
                W_reset = jnp.linalg.solve(subM.T, M.T).T
                return logabs_reset, sign_reset, W_reset
            logabs_reset, sign_reset, Ws_reset = vmap(reset_W)(xs_new, Ws_new)
            return xs_new, logabs_reset, sign_reset, Ws_reset

        def no_reset():
            return xs_new, logabs_new, sign_new, Ws_new

        state_new = jax.lax.cond(step % reset_interval == 0, do_reset, no_reset)

        return (state_new, step + 1), state_new

    keys = jax.random.split(key, n_steps)
    init_step = 0
    final_carry, hist = jax.lax.scan(single_step, (walker_state, init_step), keys)
    final_state = final_carry[0]

    return final_state, hist

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
def train(M_init, t_edges, V_edges, n_iter=100, n_walkers=64, n_steps=20, lr=0.05, \
        n_thermal=50, batch_size=1000, key=jax.random.PRNGKey(42)):
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

        # 构造当前 Slater matrix 并正交化
        M = reconstruct_M(params).astype(jnp.complex64)
        #M, _ = jnp.linalg.qr(M)
        #params = {
        #    "real": M.real,
        #    "imag": M.imag,
        #}

        # 初始化 walker
        walker_state = init_walker_states(M, N, Nf, n_walkers, key_init)

        # 热化
        walker_state, _ = sample_batch_rank1_stable(M, walker_state, n_thermal, key=key_thermal, reset_interval=50)

        # Metropolis 采样
        walker_state, (xs, _, _, Ws) = sample_batch_rank1_stable(M, walker_state, n_steps, key=key_sample, reset_interval=50)

        # 重整为 [n_samples, ...]
        xs = xs.reshape(-1, xs.shape[-1])
        Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])
        n_samples = (xs.shape[0] // batch_size) * batch_size
        xs = xs[:n_samples]
        Ws = Ws[:n_samples]

        # SR 优化
        ##params = sr_update_real(params, xs, Ws, t_edges, V_edges, lr=lr)
        params = sr_update_real_minibatch(params, xs, Ws, t_edges, V_edges, \
                batch_size=int(batch_size), lr=lr)

        # 计算能量
        flat_xs = xs.reshape(-1, xs.shape[-1])      # shape: (n_steps * n_walkers, Nf)
        flat_Ws = Ws.reshape(-1, Ws.shape[-2], Ws.shape[-1])  # shape: (n_steps * n_walkers, N, Nf)
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
    t1, t2, V1, V2 = 1.0, 0.5, 0.0, 0.0
    m = 1.00
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

    t_edges = jnp.array(t_edges)  # shape (n_t, 3)
    V_edges = jnp.array(V_edges)  # shape (n_v, 3)
    
    M_final = train(M0, t_edges, V_edges, n_iter=2, n_walkers=20, n_steps=50, \
            lr=1e-2, n_thermal=11, batch_size=1, key=key)
    return M_final

if __name__ == "__main__":
    main()
