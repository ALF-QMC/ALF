#!/usr/bin/env python3
import h5py
import numpy as np
import matplotlib.pyplot as plt

def load_density(filename):
    return np.loadtxt(filename)

def parse_v1_v2(group_name):
    parts = group_name.replace("V1_", "").split("_V2_")
    V1 = float(parts[0])
    V2 = float(parts[1])
    return V1, V2

def main():
    #target_density = load_density("density.dat")
    target_density = load_density("ndeni_scal_list")
    target_density = target_density.flatten()
    Nsite = len(target_density)

    V1_list = []
    V2_list = []
    loss_list = []

    best_diff = 1e10
    best_group = None
    best_n_mean = None

    with h5py.File("tV_HF_results.h5", "r") as f:
        for group_name in f.keys():
            group = f[group_name]
            n_mean = np.array(group["n_mean"]).reshape(-1)

            V1, V2 = parse_v1_v2(group_name)
            diff = np.mean((target_density - n_mean) ** 2)
            if (target_density[0] <= 0.50):
                n_tmp = np.array(group["n_mean"])
                n_swap = n_tmp[..., ::-1]
                n_mean = n_swap.reshape(-1)

            V1_list.append(V1)
            V2_list.append(V2)
            loss_list.append(diff)

            if diff < best_diff:
                best_diff = diff
                best_group = group
                best_group_name = group_name
                best_n_mean = n_mean

        # 打印最佳 (V1, V2)
        best_V1, best_V2 = parse_v1_v2(best_group_name)
        print(f"Best match: V1={best_V1}, V2={best_V2}, diff={best_diff:.6e}")

        # 保存最佳 slater
        slater = np.array(best_group["slater"])
        print(slater[1,:])
        #print(slater.dtype)
        slater_imag = np.zeros_like(slater)
        data = np.stack([slater.T, slater_imag.T], axis=-1)
        with h5py.File("trial_vmc.h5", "w") as f_out:
            f_out.create_dataset("phi_trial", data=data, dtype=np.float64)
        with h5py.File("M_slat.h5", "w") as f_out:
            f_out.create_dataset("M_slat", data=slater)
        print("Saved best Slater matrix to 'trial_vmc.h5'")

        # 保存 n_mean_opt.txt
        combined = np.stack([target_density, best_n_mean], axis=1)
        np.savetxt("n_mean_opt.txt", combined, fmt="%.8f", header="n_in   n_opt")
        print("Saved comparison to 'n_mean_opt.txt'")

    # 画 loss 分布
    plot_loss(V1_list, V2_list, loss_list)

def plot_loss(V1_list, V2_list, loss_list):
    V1_arr = np.array(V1_list)
    V2_arr = np.array(V2_list)
    loss_arr = np.array(loss_list)

    unique_V1 = np.sort(np.unique(V1_arr))
    unique_V2 = np.sort(np.unique(V2_arr))

    if len(unique_V1) > 1 and len(unique_V2) > 1:
        # 热点图
        V1_grid, V2_grid = np.meshgrid(unique_V1, unique_V2, indexing="xy")

        loss_grid = np.full((len(unique_V1), len(unique_V2)), np.nan, dtype=np.float64)

        for V1, V2, loss in zip(V1_arr, V2_arr, loss_arr):
            i = np.where(unique_V1 == V1)[0][0]
            j = np.where(unique_V2 == V2)[0][0]
            #loss_grid[i, j] = loss
            loss_grid[i, j] = np.log(loss)

        plt.figure(figsize=(6, 5))
        im = plt.imshow(loss_grid.T, origin="lower", extent=[unique_V1[0], unique_V1[-1], unique_V2[0], unique_V2[-1]],
                        aspect="auto")
        plt.colorbar(im, label="Loss")
        plt.xlabel("V1")
        plt.ylabel("V2")
        plt.title("Loss Heatmap")
        plt.tight_layout()
        plt.savefig("loss_heatmap.png", dpi=150)
        plt.show()
    else:
        # 曲线图
        if len(unique_V1) > 1:
            x = V1_arr
            xlabel = "V1"
        else:
            x = V2_arr
            xlabel = "V2"

        plt.figure(figsize=(6, 4))
        plt.plot(x, loss_arr, marker='o')
        plt.xlabel(xlabel)
        plt.ylabel("Loss")
        plt.title("Loss Curve")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("loss_curve.png", dpi=150)
        plt.show()

if __name__ == "__main__":
    main()
