#!/usr/bin/env python3
"""Script for automatically parsing parameters of Hamiltonian."""
# pylint: disable=invalid-name
# pylint: disable=consider-using-f-string

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2022, The ALF Project"
__license__ = "GPL"

from argparse import ArgumentParser
from pprint import pprint
import os

from parse_ham_mod import parse, create_read_write_par, get_ham_names_ham_files


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
        '--print_obj_list', action='store_true',
        help='Print list of Hamiltonian object files, derived from '
             '"Hamiltonians.list" and "Hamiltonians.list.d".'
        )
    parser.add_argument(
        '--create_hamiltonians_interface', action='store_true',
        help='Create "Hamiltonians_interface.h", derived from '
             '"Hamiltonians.list" and "Hamiltonians.list.d".'
        )
    parser.add_argument(
        '--create_hamiltonians_case', action='store_true',
        help='Create "Hamiltonians_case.h", derived from "Hamiltonians.list" '
             'and "Hamiltonians.list.d".'
        )
    parser.add_argument(
        '--create_read_write_par', action='store_true',
        help='Parse hamiltonian named in "Hamiltonians.list" and '
             '"Hamiltonians.list.d" and write '
             'Subroutines read_parameters() and write_parameters_hdf5().'
        )
    args = parser.parse_args()

    if args.test_file is not None:
        for filename in args.test_file:
            print("Parsing file: {}".format(filename))
            parameters = parse(filename)
            print("Results:")
            pprint(parameters)

    if (args.print_obj_list or args.create_hamiltonians_interface
        or args.create_hamiltonians_case or args.create_read_write_par):
        ham_names, ham_files = get_ham_names_ham_files('Hamiltonians.list')

    if args.print_obj_list:
        print(' '.join([ham_file.replace('.F90', '.o') for ham_file in ham_files]))

    if args.create_hamiltonians_interface:
        hamiltonians_interface_str = 'interface\n'
        for ham_name, ham_file in zip(ham_names, ham_files):
            hamiltonians_interface_str = f'''{hamiltonians_interface_str}
                module subroutine Ham_Alloc_{ham_name}()
                end subroutine Ham_Alloc_{ham_name}
            '''
        hamiltonians_interface_str = f'''{hamiltonians_interface_str}\nend interface'''
        with open('Hamiltonians_interface.h', 'w', encoding='UTF-8') as f:
            f.write(hamiltonians_interface_str)

    if args.create_hamiltonians_case:
        hamiltonians_case_str = ''
        for ham_name, ham_file in zip(ham_names, ham_files):
            hamiltonians_case_str = f'''{hamiltonians_case_str}
                Case ('{ham_name.upper()}')
                    call Ham_Alloc_{ham_name}'''
        with open('Hamiltonians_case.h', 'w', encoding='UTF-8') as f:
            f.write(hamiltonians_case_str)

    if args.create_read_write_par:
        for ham_name, ham_file in zip(ham_names, ham_files):
            parameters = parse(ham_file)
            # pprint(parameters)
            filename = os.path.join(
                os.path.dirname(ham_file),
                'Hamiltonian_{}_read_write_parameters.F90'.format(ham_name)
                )
            print('filenames:', ham_file, filename)
            create_read_write_par(filename, parameters, ham_name)
