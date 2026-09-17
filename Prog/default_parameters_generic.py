"""Generic (Hamiltonian independent) parameters with default values."""

__author__ = "Jonas Schwab"
__copyright__ = "Copyright 2022, The ALF Project"
__license__ = "GPL"

from collections import OrderedDict


# Generic (Hamiltonian independent) ALF parameters with default values.
_PARAMS_GENERIC = OrderedDict([
    ('VAR_QMC',
        {'CPU_MAX': {'comment': 'Code stops after CPU_MAX hours, if 0 or '
                                'not specified, the code stops after '
                                'Nbin bins',
                     'value': 0.0},
         'Delta_t_Langevin_HMC': {'comment': 'Time step for Langevin or '
                                             'HMC',
                                  'value': 0.1},
         'Global_moves': {'comment': 'Allows for global moves in space '
                                     'and time.',
                          'value': False},
         'Global_tau_moves': {'comment': 'Allows for global moves on a '
                                         'single time slice.',
                              'value': False},
         'HMC': {'comment': 'HMC update', 'value': False},
         'LOBS_EN': {'comment': 'End measurements at time slice LOBS_EN',
                     'value': 0},
         'LOBS_ST': {'comment': 'Start measurements at time slice '
                                'LOBS_ST',
                     'value': 0},
         'sequential': {'comment': 'Conventional updating scheme', 'value': True},
         'Langevin': {'comment': 'Langevin update', 'value': False},
         'Leapfrog_steps': {'comment': 'Number of leapfrog iterations',
                            'value': 0},
         'N_HMC_sweeps': {'comment': 'Number of HMC sweeps',
                            'value': 1},
         'Amplitude':    {'comment': 'Width of the box distribution for update of '
                                     'type t=3,4 fields.  Defaults  to 1.0.',
                            'value': 1.0},
         'Ltau': {'comment': '1 to calculate time-displaced Green '
                             'functions; 0 otherwise.',
                  'value': 1},
         'Max_Force': {'comment': 'Max Force for Langevin', 'value': 5.0},
         'Sequential_MALA': {'comment': 'For MALA updates', 'value': False},
         'Max_Force_MALA_sequential': {'comment': 'Max Force for '
                                                  'sequential MALA '
                                                  'updates',
                                      'value': 1.0},
         'Delta_t_MALA_sequential': {'comment': 'Time step for '
                                                'sequential MALA '
                                                'updates',
                                    'value': 0.1},
         'Global_tau_MALA_moves': {'comment': 'For Global_tau_MALA_move',
                                   'value': False},
         'Global_MALA_moves': {'comment': 'For Global (space-time) MALA updates',
                                   'value': False},
         'N_Global_tau_MALA': {'comment': 'Number of global tau MALA '
                                          'moves that will be carried '
                                          'out on a single time slice.',
                               'value': 1},
         'Delta_t_MALA_global_tau': {'comment': 'Time step for global '
                                                'tau MALA updates',
                                     'value': 0.1},
         'Max_Force_MALA_global_tau': {'comment': 'Max Force for global '
                                                  'tau MALA updates',
                                       'value': 1.0},
         'N_global': {'comment': 'Number of global moves per sweep.',
                      'value': 1},
         'N_global_tau': {'comment': 'Number of global moves that will '
                                     'be carried out on a single time '
                                     'slice.',
                          'value': 1},
         'Nbin': {'comment': 'Number of bins.', 'value': 5},
         'Nsweep': {'comment': 'Number of sweeps per bin.', 'value': 20},
         'Nt_sequential_end': {'comment': '', 'value': -1},
         'Nt_sequential_start': {'comment': '', 'value': 0},
         'Nwrap': {'comment': 'Stabilization. Green functions will be '
                              'computed from scratch after each time '
                              'interval Nwrap*Dtau.',
                   'value': 10},
         'Propose_S0': {'comment': 'Proposes single spin flip moves with '
                                   'probability exp(-S0).',
                        'value': False}}),
    ('VAR_errors',
        {'N_Back': {'comment': 'If set to 1, substract background in '
                               'correlation functions. Is ignored in '
                               'Python analysis.',
                    'value': 1},
         'N_Cov': {'comment': 'If set to 1, covariance computed for '
                              'time-displaced correlation functions. '
                              'Is ignored in Python analysis.',
                   'value': 0},
         'N_auto': {'comment': 'If > 0, calculate autocorrelation. '
                                'Is ignored in Python analysis.',
                    'value': 0},
         'N_rebin': {'comment': 'Rebinning: Number of bins to combine '
                                'into one.',
                     'value': 1},
         'N_skip': {'comment': 'Number of bins to be skipped.',
                    'value': 1},
         'N_BZ_Zones': {'comment': 'Number of Brillouin zones that will '
                        'be covered if extended zone scheme is activated.',
                        'value': 1},
         'Extended_Zone': {'comment': 'If true carries out the fourier '
                           'transform in the extened zone scheme.',
                           'value': False},
        }),
    ('VAR_TEMP',
        {'N_Tempering_frequency': {'comment': 'The frequency, in units '
                                              'of sweeps, at which the '
                                              'exchange moves are '
                                              'carried out.',
                                   'value': 10},
         'N_exchange_steps': {'comment': 'Number of exchange moves.',
                              'value': 6},
         'Tempering_calc_det': {'comment': 'Specifies whether the '
                                           'fermion weight has to be '
                                           'taken into account while '
                                           'tempering. Can be set to '
                                           '.F. if the parameters '
                                           'that get varied only '
                                           'enter the Ising action S_0',
                                'value': True},
         'mpi_per_parameter_set': {'comment': 'Number of mpi-processes '
                                              'per parameter set.',
                                   'value': 2}}),
    ('VAR_Max_Stoch',
        {'Checkpoint': {'comment': '', 'value': False},
         'NBins': {'comment': 'Number of bins for Monte Carlo.',
                   'value': 250},
         'NSweeps': {'comment': 'Number of sweeps per bin.', 'value': 70},
         'N_alpha': {'comment': 'Number of temperatures.', 'value': 14},
         'Ndis': {'comment': 'Number of boxes for histogram.',
                  'value': 2000},
         'Ngamma': {'comment': 'Number of Dirac delta-functions for '
                               'parametrization.',
                    'value': 400},
         'Nwarm': {'comment': 'The Nwarm first bins will be omitted.',
                   'value': 20},
         'Om_en': {'comment': 'Frequency range upper bound.',
                   'value': 10.0},
         'Om_st': {'comment': 'Frequency range lower bound.',
                   'value': -10.0},
         'R': {'comment': '', 'value': 1.2},
         'Stochastic': {'comment': 'If true, then stochastic MaxEnt is '
                                   'used. If false, then classic '
                                   'MaxEnt is used.',
                        'value': True},
         'Tolerance': {'comment': '', 'value': 0.1},
         'alpha_st': {'comment': '', 'value': 1.0}})
    ])
