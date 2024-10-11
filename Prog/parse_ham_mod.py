"""
Collection of functions for automatically parsing parameters of Hamiltonian.

Includes functions for writing subroutines read_parameters() and
write_parameters_hdf5().
"""
# pylint: disable=invalid-name
# pylint: disable=consider-using-f-string

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2022, The ALF Project"
__license__ = "GPL"

import glob


def get_ham_names_ham_files(ham_list_file):
    """Return list of Hamiltonian names and files.
    
    Also parses files in directory `ham_list_file` + '.d',
    so usually ${ALF_DIR}/Prog/Hamiltonians.list and files in
    ${ALF_DIR}/Prog/Hamiltonians.list.d
    """
    ham_list_files = [ham_list_file, *glob.glob(ham_list_file + '.d/*')]

    lines = []
    for file in ham_list_files:
        with open(file, 'r', encoding='UTF-8') as f:
            for line in f.read().splitlines():
                if not line.strip()[0] == '#':
                    lines.append(line.strip().split())

    ham_names = []
    ham_files = []
    for line in lines:
        ham_names.append(line[0])
        try:
            ham_files.append(line[1])
        except IndexError:
            ham_files.append(
                'Hamiltonians/Hamiltonian_{}_smod.F90'.format(line[0]))
    return ham_names, ham_files

def parse(filename):
    """Parse Fortran file for parameter lists.

    Each list has to start with a line containing '#PARAMETERS START#' and end
    with a line containing '#PARAMETERS END#'. The name written in the same
    line as '#PARAMETERS START#' is the namelist name.
    """
    parameters = {}
    with open(filename, 'r', encoding='UTF-8') as f:
        lines = f.readlines()

    do_parse = False
    for i, line in enumerate(lines):
        if '#PARAMETERS START#' in line:
            do_parse = True
            name_key = line.split('#PARAMETERS START#')[1].strip()
            if not name_key.isidentifier():
                raise Exception(
                    'In line {} of {}: '.format(i+1, filename) +
                    '"{}" is not a valid namelist name.'.format(name_key)
                    )
            namelist = {}
            parameters[name_key] = namelist
            continue
        if '#PARAMETERS END#' in line:
            do_parse = False
            continue
        if do_parse and '::' in line:
            par = parse_line(line)
            par_name = par.pop('name')
            namelist[par_name] = par
    return parameters


def parse_line(line):
    """Parse single line in Fortran file for parameter."""
    parameter = {}
    dtype, rest = line.split('::', maxsplit=1)
    dtype = dtype.strip().lower()
    if '!' in rest:
        assignment, comment = rest.strip().split('!', maxsplit=1)
        assignment = assignment.strip()
        parameter['comment'] = comment.strip(' !')
    else:
        assignment = rest.strip()
        parameter['comment'] = ''
    name, value = assignment.split('=')
    # parameter['name'] = name.split('(')[0].strip()
    parameter['name'] = name.strip()
    value = value.strip()
    parameter['defined_in_base'] = dtype.startswith('!')
    if '[' in value:
        parameter['value'] = []
        for val in value.strip('[]').split(','):
            parameter['value'].append(_to_value(dtype, val))
    else:
        parameter['value'] = _to_value(dtype, value)

    return parameter


def _max_len(dictionary):
    return max(len(i) for i in dictionary)


def _dtype_name(parameter):
    if isinstance(parameter, list):
        return _dtype_name(parameter[0])
    if isinstance(parameter, bool):
        return 'logical'
    if isinstance(parameter, float):
        return 'real(dp)'
    # if isinstance(parameter, complex):
    #     return 'complex(dp)'
    if isinstance(parameter, int):
        return 'integer'
    if isinstance(parameter, str):
        return 'Character(len=64)'
    raise Exception('Error in "_dtype_name": unrecognized type')


def _convert_par_to_str(parameter):
    """Converts a given parameter value to a string that can be
    written into a parameter file.
    """
    if isinstance(parameter, list):
        pars = [_convert_par_to_str(p) for p in parameter]
        return '[{}]'.format(', '.join(pars))
    if isinstance(parameter, bool):
        if parameter:
            return '.true.'
        return '.false.'
    if isinstance(parameter, float):
        if 'e' in '{}'.format(parameter):
            return '{}'.format(parameter).replace('e', 'd')
        return '{}d0'.format(parameter)
    # if isinstance(parameter, complex):
    #     return '({}, {})'.format(_convert_par_to_str(parameter.real),
    #                              _convert_par_to_str(parameter.imag))
    if isinstance(parameter, int):
        return '{}'.format(parameter)
    if isinstance(parameter, str):
        return '"{}"'.format(parameter)

    raise Exception('Error in "_convert_par_to_str": unrecognized type')


def create_read_write_par(filename, parameters, ham_name):
    """Write subroutines read_parameters() and write_parameters_hdf5()."""
    INDENT = 9     # Number of indentation spaces
    LINE_MAX = 70  # Maximal line length
    TEMPLATE = """! This file is automatically generated by parse_ham.py
! Please note that all manual changes will be overwritten during compilation!
! Supplies subroutine for automatically reading and writing the parameters of "##HAM_NAME##"

      Subroutine read_parameters()
         use iso_fortran_env, only: error_unit
#ifdef MPI
         Use mpi
#endif
         Implicit none

         !Local variables
         integer            :: ierr, unit_para
         Character (len=64) :: file_info, file_para
         Integer, parameter :: dp = kind(0.d0)

         ##NAMELIST##

#ifdef MPI
         Integer        :: Isize, Irank, igroup, irank_g, isize_g
         Integer        :: STATUS(MPI_STATUS_SIZE)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         IF (ierr /= 0) THEN
            WRITE(error_unit,*) 'read_parameters: Error in MPI', ierr
            Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
         END IF
         igroup           = irank/isize_g
#endif

         ##PARAMETER_DEF##
#ifdef MPI
         If (Irank_g == 0 ) then
#endif
#ifdef TEMPERING
            write(file_para,'(A,I0,A)') "Temp_", igroup, "/parameters"
            write(file_info,'(A,I0,A)') "Temp_", igroup, "/info"
#else
            file_para = "parameters"
            file_info = "info"
#endif

            !Read in parameters
            OPEN(NEWUNIT=unit_para, FILE=file_para, STATUS='old', ACTION='read', IOSTAT=ierr)
            IF (ierr /= 0) THEN
               WRITE(error_unit,*) 'unable to open <parameters>', ierr
               Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
            END IF
            ##READ_VAR##
            close(unit_para)

#ifdef MPI
         Endif
         !Broadcast parameters to all MPI tasks
         ##MPI_BCAST##
#endif

      end subroutine read_parameters


#ifdef HDF5
      Subroutine write_parameters_hdf5(filename)
         Use hdf5
         use alf_hdf5
         Implicit none

         Character (len=64), intent(in) :: filename

         Logical            :: link_exists
         INTEGER(HID_T)     :: file_id, group_id, group_id2
         INTEGER            :: ierr

         ! Open HDF5 file
         CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, ierr)

         CALL h5lexists_f(file_id, "parameters", link_exists, ierr)
         if ( .not. link_exists ) then
           ! Create parameters-group
           call h5gcreate_f(file_id, "parameters", group_id, ierr)
         else
           ! Open parameters-group
           call h5gopen_f (file_id, "parameters", group_id, ierr)
         endif

         ##TEST_ATTRS##

         call h5gclose_f(group_id, ierr)
         call h5fclose_f(file_id, ierr)

      end Subroutine write_parameters_hdf5
#endif
"""

    f = open(filename, 'w', encoding='UTF-8')
    for line in TEMPLATE.splitlines(keepends=True):
        if '##NAMELIST##' in line:
            for nlist_name, nlist in parameters.items():
                s = '{}NAMELIST /{}/  '.format(INDENT*' ', nlist_name)
                for par_name in nlist:
                    par_name = par_name.split('(')[0]
                    if len(s) < LINE_MAX:
                        s = '{}{}, '.format(s, par_name)
                    else:
                        f.write('{}&\n'.format(s))
                        s = '{}     &     {}, '.format(INDENT*' ', par_name)
                f.write(s[:-2]+'\n')
        elif '##PARAMETER_DEF##' in line:
            for nlist_name, nlist in parameters.items():
                f.write('{}!Parameters {}\n'.format(INDENT*' ', nlist_name))

                names_len = _max_len(nlist)
                # dtypes_str = [_dtype_name(par['value']) for
                #               par_name, par in nlist.items()]
                # dtypes_len = _max_len(dtypes_str)
                pars_str = [_convert_par_to_str(par['value']) for
                            par_name, par in nlist.items()]
                pars_len = _max_len(pars_str)
                # comments = [par[1] for par_name, par in nlist.items()]
                # for i in range(len(nlist)):

                for par_name, par in nlist.items():
                    # s = '{}{} :: {} = {} !'.format(
                    #     INDENT*' ',
                    #     _dtype_name(par['value']).ljust(dtypes_len),
                    #     par_name.ljust(names_len),
                    #     _convert_par_to_str(par['value']).ljust(pars_len)
                    #     )
                    s = '{}{} = {} !'.format(
                        INDENT*' ',
                        par_name.ljust(names_len),
                        _convert_par_to_str(par['value']).ljust(pars_len)
                        )
                    comment = par['comment'].split(' ')
                    comment_indent = min(40, len(s)-2)
                    for word in comment:
                        if len(s) < LINE_MAX:
                            s = '{} {}'.format(s, word)
                        else:
                            f.write('{} \n'.format(s))
                            s = '{} ! {}'.format(comment_indent*' ', word)
                    f.write(s)
                    f.write('\n')
        elif '##READ_VAR##' in line:
            for nlist_name, nlist in parameters.items():
                f.write('{}   REWIND(unit_para)\n'.format(INDENT*' '))
                s = '{}   READ(unit_para, NML={})\n'.format(
                    INDENT*' ', nlist_name)
                f.write(s)

        elif '##MPI_BCAST##' in line:
            fstring = '{}CALL MPI_BCAST({},{:>3},{},0,Group_Comm,ierr)\n'
            for nlist_name, nlist in parameters.items():
                names_len = _max_len(nlist)
                for par_name, par in nlist.items():
                    s = fstring.format(
                        INDENT*' ',
                        par_name.split('(')[0].ljust(names_len),
                        _get_mpi_len(par['value']),
                        _get_mpi_dtype(par['value']),
                        )
                    f.write(s)
        elif '##TEST_ATTRS##' in line:
            for nlist_name, nlist in parameters.items():
                s = """CALL h5lexists_f(group_id, "{}", link_exists, ierr)
                       if ( .not. link_exists ) then
                         call h5gcreate_f(group_id, "{}", group_id2, ierr)
                         call h5gclose_f(group_id2, ierr)
                       endif\n""".format(nlist_name.lower(), nlist_name.lower())
                f.write(s)
                for par_name, par in nlist.items():
                    s = '{}call test_attribute(group_id, "{}", "{}", {}, ierr)\n'.format(
                        INDENT*' ', nlist_name.lower(), par_name.lower(), par_name
                        )
                    f.write(s)
        else:
            f.write(line.replace('##HAM_NAME##', ham_name))
    f.close()


def _to_value(dtype, value):
    if 'real' in dtype:
        return float(value.replace('d', 'e').replace('D', 'e'))
    # if 'complex' in dtype:
    #     tmp = value.strip('()').split(',')
    #     if not len(tmp) == 2:
    #         raise Exception(
    #             '"{}", "{}" cannot be mapped to complex'.format(dtype, value))
    #     return complex(float(tmp[0]), float(tmp[1]))
    if 'integer' in dtype:
        return int(value)
    if 'character' in dtype:
        return value.strip('"\'')
    if 'logical' in dtype:
        if 't' in value or 'T' in value:
            return True
        if 'f' in value or 'F' in value:
            return False

        raise Exception(
            '"{}" can not be mapped to bool.'.format(value))

    raise Exception('"{}" can not be mapped to a type'.format(dtype))


def _get_mpi_dtype(parameter):
    if isinstance(parameter, list):
        return _get_mpi_dtype(parameter[0])
    if isinstance(parameter, bool):
        return 'MPI_LOGICAL  '
    if isinstance(parameter, float):
        return 'MPI_REAL8    '
    # if isinstance(parameter, complex):
    #     return 'MPI_COMPLEX16'
    if isinstance(parameter, int):
        return 'MPI_INTEGER  '
    if isinstance(parameter, str):
        return 'MPI_CHARACTER'
    raise Exception('Error in "_get_mpi_dtype": unrecognized type')


def _get_mpi_len(parameter):
    if isinstance(parameter, list):
        return _get_mpi_len(parameter[0]) * len(parameter)
    if isinstance(parameter, (bool, float, complex, int)):
        return 1
    if isinstance(parameter, str):
        return 64
    raise Exception('Error in "_get_mpi_len": unrecognized type')
