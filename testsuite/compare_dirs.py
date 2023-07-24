#!/usr/bin/env python3
"""Testing two directories of ALF simulationsfor identical results."""
# pylint: disable=invalid-name

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2023, The ALF Project"
__license__ = "GPL"

import sys
import argparse

import numpy as np
from py_alf.ana import load_res
from py_alf.analysis import analysis
from py_alf.utils import find_sim_dirs

def analyze_and_load(root_dir):
    dirs = find_sim_dirs(root_dir)
    for dir in dirs:
        analysis(dir)
    return load_res(dirs)


def compare_dirs(dir_R, dir_T, resultsfile):
    obs_R = analyze_and_load(dir_R)
    obs_T = analyze_and_load(dir_T)

    with open(resultsfile, 'w', encoding='UTF-8') as f:
        test_all = True
        for name in obs_R:
            test = True
            for dat_R, dat_T in zip(obs_R[name], obs_T[name]):
                try:
                    test_temp = np.allclose(dat_R, dat_T)
                except TypeError:
                    pass
                test = test and test_temp
            f.write(f'{name}: {test}\n')
            test_all = test_all and test
        f.write(f'\nTotal: {test}\n')
    return test_all


def _get_arg_parser():
    parser = argparse.ArgumentParser(
        description='Testing two directories of ALF simulations'
                    'for identical results.')
    parser.add_argument(
        'dir_R', help='Directory containing reference simulation.')
    parser.add_argument(
        'dir_T', help='Directory containing test simulation.')
    parser.add_argument(
        '--results', default="test.txt",
        help='Filename to store results in.   (default: "test.txt")')
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()

    test_success = compare_dirs(args.dir_R, args.dir_T, args.results)

    if test_success:
        print("Test sucessful")
        sys.exit(0)
    else:
        print("Test failed")
        with open(args.results, 'r', encoding='UTF-8') as f:
            print(f.read())
        sys.exit(1)
