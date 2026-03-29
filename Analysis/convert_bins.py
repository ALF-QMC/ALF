#!/usr/bin/env python3
"""Copies parameters from parameters file to hdf5 file.
"""
import sys
import os
import subprocess
from argparse import ArgumentParser

import h5py
import f90nml

def get_val(default_parameters, nml, nlist_name, par_name):
    r"""Get value from namelist or return default value."""
    nlist_name = nlist_name.lower()
    par_name = par_name.lower()
    try:
        val = nml[nlist_name][par_name]
    except KeyError:
        val = default_parameters[nlist_name][par_name]['value']
    return val


def copy_parameters(sim_dir, hamiltonian_file):
    r"""Copy parameters from parameters file and hamiltonian to HDF5 file."""
    nml = f90nml.read(os.path.join(sim_dir, 'parameters'))

    default_parameters = parse(hamiltonian_file)
    filename = os.path.join(sim_dir, 'data.h5')

    # Fix for Mz=true
    if (get_val(default_parameters, nml, 'var_ham_name', 'ham_name').lower() == 'hubbard' and
       get_val(default_parameters, nml, 'var_hubbard', 'mz')):
        print('Applying Mz=true fix')
        mz_fix = True
    else:
        mz_fix = False

    with h5py.File(filename, 'w') as f:
        for nlist_name, nlist in default_parameters.items():
            nlist_name = nlist_name.lower()
            groupname = f"parameters/{nlist_name}"
            f.create_group(groupname)
            for par_name in nlist.keys():
                par_name = par_name.lower()
                val = get_val(default_parameters, nml, nlist_name, par_name)
                if mz_fix and nlist_name == 'var_model_generic':
                    # Fix for Mz=true
                    if par_name == 'n_sun':
                        val = val // 2
                    if par_name == 'n_fl':
                        val = 2

                if isinstance(val, bool):
                    val = int(val)
                
                if isinstance(val, str):
                    f[groupname].attrs.create(
                        par_name, val,
                        dtype=h5py.string_dtype(encoding='ascii', length=len(val)))
                else:
                    f[groupname].attrs.create(par_name, val)


def _get_arg_parser():
    parser = ArgumentParser(
        description='Convert plain text bins and parameters to HDF5 file. '
        'Moves plain text bin and associated _info files to "old_bins" '
        'subfolder or optionally removes them. '
        'One could run a bigger batch with '
        r'`find . -name "Ener_scal" -execdir ${ALF_DIR}/Analysis/copy_parameters.py \;` '
        'Where "Ener_scal" can be replaced by any file name for recognizing simulation folders.',
        )
    parser.add_argument(
        '--alfdir', default=os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
        help="Path to ALF directory. Otherwise, it is determined by location of this script.")
    parser.add_argument(
        '--remove-old-bins', action='store_true',
        help='Remove old plain text bins instead of moving them to "old_bins" subfolder.')
    parser.add_argument(
        'directories', nargs='*',
        help='Directories')
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()

    directories = args.directories if args.directories else ['.']

    sys.path.append(os.path.join(args.alfdir, 'Prog'))
    from parse_ham_mod import parse, get_ham_names_ham_files

    convert_scal = os.path.join(args.alfdir, 'Analysis', 'convert_scal.out')
    convert_latt = os.path.join(args.alfdir, 'Analysis', 'convert_latt.out')
    convert_local = os.path.join(args.alfdir, 'Analysis', 'convert_local.out')

    ham_names, ham_files = get_ham_names_ham_files(
        os.path.join(args.alfdir, 'Prog', 'Hamiltonians.list'))
    ham_file_dict = {}
    for ham_name, ham_file in zip(ham_names, ham_files):
        ham_file_dict[ham_name] = ham_file

    for directory in directories:
        nml = f90nml.read(os.path.join(directory, 'parameters'))
        ham_name = nml['VAR_ham_name']['ham_name']
        ham_file = os.path.join(args.alfdir, 'Prog', ham_file_dict[ham_name])
        copy_parameters(directory, ham_file)
        print(f'Convert bins in folder "{os.path.abspath(directory)}"')
        os.makedirs(os.path.join(directory, 'old_bins'), exist_ok=True)
        for dir_content in os.listdir(path=directory):
            old_bins = []
            if dir_content.endswith('_scal'):
                subprocess.run([convert_scal, dir_content, 'data.h5'], check=True, cwd=directory)
                old_bins.append(dir_content)
            elif dir_content.endswith('_eq') or dir_content.endswith('_tau'):
                subprocess.run([convert_latt, dir_content, 'data.h5'], check=True, cwd=directory)
                old_bins.append(dir_content)
            elif dir_content.endswith('_local') or dir_content.endswith('_localtau'):
                subprocess.run([convert_local, dir_content, 'data.h5'], check=True, cwd=directory)
                old_bins.append(dir_content)

            info_files = [f+'_info' for f in old_bins if os.path.isfile(os.path.join(directory, f+'_info'))]
            for file in old_bins+info_files:
                if args.remove_old_bins:
                    os.remove(os.path.join(directory, file))
                else:
                    os.rename(os.path.join(directory, file),
                        os.path.join(directory, 'old_bins', file))
