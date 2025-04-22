#!/usr/bin/env python3
import jax
import jax.numpy as jnp
import h5py
import numpy as np
import csv

jax.config.update("jax_enable_x64", True)


def compute_G(M):
    Q, _ = jnp.linalg.qr(M)
    return Q @ Q.conj().T

def main(den_file="den.dat", phi_h5="all_phi_trials.h5", info_csv="results_summary.csv"):
    # 读取目标密度分布 n_in
    n_in = np.loadtxt(den_file)
    n_in = jnp.array(n_in).flatten()

    best_loss = float("inf")
    best_v1 = None
    best_n = None

    # 读取 CSV 获取所有 V1 和对应键
    entries = []
    with open(info_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            entries.append((float(row['V1']), row['filename']))

    # 打开统一的 HDF5 文件
    with h5py.File(phi_h5, "r") as f:
        for v1, key in entries:
            try:
                data = f[key][:]  # 读取 (ndim, n_part, 2)
                M_real = data[..., 0]
                M_imag = data[..., 1]
                M = jnp.array(M_real + 1j * M_imag)
                G = compute_G(M)
                n = jnp.real(jnp.diag(G))
                loss = jnp.mean((n - n_in) ** 2)
                print(f"V1={v1:.3f}, loss={loss:.6e}")
                if loss < best_loss:
                    best_loss = float(loss)
                    best_v1 = v1
                    best_n = n
            except Exception as e:
                print(f"Warning: failed to process {key}: {e}")

    print("\nBest match:")
    print(f"  V1 = {best_v1:.6f}")
    print(f"  Loss = {best_loss:.6e}")
    np.savetxt("best_density_vs_input.txt", jnp.stack([n_in, best_n], axis=1), header="n_in  n_pred")


if __name__ == '__main__':
    main()

