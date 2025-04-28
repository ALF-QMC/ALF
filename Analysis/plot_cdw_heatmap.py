#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import re

def parse_v1_v2(key):
    match = re.match(r"V1_(\d+\.\d+)_V2_(\d+\.\d+)", key)
    if match:
        v1 = float(match.group(1))
        v2 = float(match.group(2))
        return v1, v2
    else:
        return None, None

def load_cdw_data(h5file_path):
    v1_list, v2_list, cdw_list = [], [], []
    with h5py.File(h5file_path, "r") as f:
        for key in f.keys():
            v1, v2 = parse_v1_v2(key)
            if v1 is not None:
                cdw = f[key]["stagger_cdw"][()]
                v1_list.append(v1)
                v2_list.append(v2)
                cdw_list.append(cdw)
    return np.array(v1_list), np.array(v2_list), np.array(cdw_list)

def plot_cdw_vs_single_param(param, cdw, param_name):
    sorted_idx = np.argsort(param)
    param_sorted = param[sorted_idx]
    cdw_sorted = cdw[sorted_idx]

    plt.figure(figsize=(8,5))
    plt.plot(param_sorted, cdw_sorted, marker='o')
    plt.xlabel(param_name)
    plt.ylabel("Staggered CDW Order")
    plt.title(f"CDW Order vs {param_name}")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_cdw_heatmap(v1, v2, cdw):
    v1_unique = np.unique(v1)
    v2_unique = np.unique(v2)
    V1_grid, V2_grid = np.meshgrid(v1_unique, v2_unique, indexing='ij')
    CDW_grid = np.full(V1_grid.shape, np.nan)

    for i in range(len(v1)):
        idx_v1 = np.where(v1_unique == v1[i])[0][0]
        idx_v2 = np.where(v2_unique == v2[i])[0][0]
        CDW_grid[idx_v1, idx_v2] = cdw[i]

    plt.figure(figsize=(8,6))
    cmap = plt.get_cmap('viridis')
    c = plt.pcolormesh(V1_grid, V2_grid, CDW_grid, shading='auto', cmap=cmap)
    plt.colorbar(c, label='Staggered CDW Order')
    plt.xlabel('V1')
    plt.ylabel('V2')
    plt.title('CDW Order Parameter Heatmap')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    h5file_path = "tV_HF_results.h5"
    v1, v2, cdw = load_cdw_data(h5file_path)

    if len(v1) == 0:
        print("No data found in HDF5 file.")
    elif len(np.unique(v2)) == 1:
        plot_cdw_vs_single_param(v1, cdw, "V1")
    elif len(np.unique(v1)) == 1:
        plot_cdw_vs_single_param(v2, cdw, "V2")
    else:
        plot_cdw_heatmap(v1, v2, cdw)
