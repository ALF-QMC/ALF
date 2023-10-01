#!/usr/bin/env python3
import sys

import yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from py_alf import analysis
from py_alf.utils import find_sim_dirs
from py_alf.ana import load_res

with open("spec.yaml", 'r', encoding='UTF-8') as f:
    spec = yaml.load(f, yaml.Loader)

dirs = find_sim_dirs()
for directory in dirs:
    analysis(directory)
res = load_res(dirs)
res.to_json("results.json")

def func(x, y0, a):
    return y0 + a*x**2
x = res.dtau
y = res.Ener_scal0
dy = res.Ener_scal0_err
popt, pcov = curve_fit(func, x, y, sigma=dy, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))

plt.figure()
plt.title(f"Test: {spec['test_name']}; Environment: {spec['env_name']}")
plt.xlabel(r'$\Delta\tau^2$')
plt.ylabel('Energy')
plt.errorbar(x**2, y, dy, fmt='.', label="Monte carlo results")
plt.plot(0, spec['ed_energy'], 'x', label='ED result')
xs = np.linspace(0., x.max())
p = plt.plot(xs**2, func(xs, *popt))
plt.errorbar(0, popt[0], perr[0], label="Extrapolated value", color=p[0].get_color())
plt.legend()
plt.savefig('test.png')

delta = popt[0] - spec['ed_energy']
sigma = abs(delta/perr[0])

print(f"Reference value: {spec['ed_energy']}")
print(f"Extrapolated value: {popt[0]}+-{perr[0]}")
print(f"Absolute deviation: {delta}")
print(f"Deviation/error: {sigma}")
test_failed = False
if abs(delta) > spec["max_delta"]:
    print("Deviation too big, test failed!")
    test_failed = True
if sigma > spec["max_sigma"]:
    print("Deviation/error too big, test failed!")
    test_failed = True
if test_failed:
    sys.exit(1)
