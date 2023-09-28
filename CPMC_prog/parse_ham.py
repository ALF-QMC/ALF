#!/usr/bin/env python3
"""Script for automatically parsing parameters of Hamiltonian."""
# pylint: disable=invalid-name
# pylint: disable=consider-using-f-string

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2022, The ALF Project"
__license__ = "GPL"

from argparse import ArgumentParser
from pprint import pprint

from parse_ham_mod import parse, create_read_write_par


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Script for parsing parameters of Hamiltonian.',
        )
    parser.add_argument(
        '--test_file', nargs='+', default=None,
        help='Test parameter specifications within this file(s). '
             'Result get printed for manual check.'
        )
    parser.add_argument(
        '--create_read_write_par', action='store_true',
        help='Parse hamiltonian named in "Hamiltonians.list" and write '
             'Subroutines read_parameters() and write_parameters_hdf5().'
        )
    args = parser.parse_args()

    if args.test_file is not None:
        for filename in args.test_file:
            print("Parsing file: {}".format(filename))
            parameters = parse(filename)
            print("Results:")
            pprint(parameters)

    if args.create_read_write_par:
        with open('Hamiltonians.list', 'r', encoding='UTF-8') as f:
            ham_names = f.read().splitlines()

        for ham_name in ham_names:
            filename = 'Hamiltonians/Hamiltonian_{}_smod.F90'.format(ham_name)
            print('filename:', filename)

            parameters = parse(filename)
            # pprint(parameters)
            filename = 'Hamiltonians/' + \
                'Hamiltonian_{}_read_write_parameters.F90'.format(ham_name)
            create_read_write_par(filename, parameters, ham_name)
