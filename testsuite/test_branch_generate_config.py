#!/usr/bin/env python3
import sys
import json
import copy
import yaml

pipeline_config = yaml.load("""
default:
  tags:
    - k8s
  artifacts:
    expire_in: 1 day

stages:
  - compile
  - simulation
#  - summarize
#
#summarize:
#  image: git.physik.uni-wuerzburg.de:25812/z03/pdi/debian:bookworm-buildd
#  stage: summarize
#  when: always
#  script:
#    - apt-get update && apt-get install -y python3 
#      #python3-yaml python3-numpy python3-pandas python3-matplotlib python3-scipy
#    - 
#  artifacts:
#    paths:
#      - ALF_data
#    when: always
#    expire_in: 1 week

.compile_template:
  stage: compile
  script:
    - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
    - export PATH="$HOME/.local/bin:$PATH"
    - git remote set-branches origin master $CI_COMMIT_BRANCH
    - git fetch --depth=1
    - git checkout $CI_COMMIT_BRANCH
    - git checkout .
    - pip install --no-deps pyALF
    - export ALF_DIR="$PWD"
    - . ./configure.sh $MACHINE noMPI HDF5 NO-INTERACTIVE
    - alf_test_branch
      --sim_pars $ALF_DIR/testsuite/test_branch_parameters.json
      --machine $MACHINE
      --mpi
      --branch_T $CI_COMMIT_BRANCH
      --branch_R master
      --devel
      --no_sim
      --no_analyze
    - rdfind -minsize 100000 -makesymlinks true ALF_data
    - symlinks -c ALF_data
  artifacts:
    paths:
      - ALF_data

.simulation_template:
  stage: simulation
  script:
    - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
    - START_DIR=$PWD
    - export PATH="$HOME/.local/bin:$PATH"
    - pip install --no-deps pyALF
    - cd ${START_DIR}/ALF_data/${TEST_NAME}
    - mpiexec -n 4 ./ALF.out
    - cd ${START_DIR}/ALF_data/${TEST_NAME}_test
    - mpiexec -n 4 ./ALF.out
    - ${START_DIR}/testsuite/compare_dirs.py
      --results ${START_DIR}/ALF_data/${TEST_NAME}.txt
      ${START_DIR}/ALF_data/${TEST_NAME}
      ${START_DIR}/ALF_data/${TEST_NAME}_test
  artifacts:
    paths:
      - ALF_data/${TEST_NAME}
      - ALF_data/${TEST_NAME}_test
      - ALF_data/${TEST_NAME}.txt
    exclude:
      - ALF_data/${TEST_NAME}/**/confout_*.h5
      - ALF_data/${TEST_NAME}/**/ALF.out
      - ALF_data/${TEST_NAME}_test/**/confout_*.h5
      - ALF_data/${TEST_NAME}_test/**/ALF.out
    when: always
    expire_in: 1 week

.analysis_template:
  stage: analyze
  script:
    - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
    - export PATH="$HOME/.local/bin:$PATH"
    - git remote set-branches origin master $CI_COMMIT_BRANCH
    - git fetch --depth=1
    - pip install --no-deps pyALF
    - export ALF_DIR="$PWD"
    - alf_test_branch
      --sim_pars $ALF_DIR/testsuite/test_branch_parameters.json
      --machine $MACHINE
      --mpi
      --branch_T $CI_COMMIT_BRANCH
      --branch_R master
      --devel
      --no_prep
      --no_sim
  artifacts:
    paths:
      - ALF_data
      - test.txt
    exclude:
      - ALF_data/**/res
    when: always
    expire_in: 1 week
""", yaml.Loader)


ENVIRONMENTS = yaml.load("""
Bullseye:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye
    variables: {MACHINE: GNU}
Bookworm:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm
    variables: {MACHINE: GNU}
#Intel21:
#    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-intel
#    variables: {MACHINE: INTEL}
Intel-2024.2:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm-intel-2024.2
    variables: {MACHINE: INTEL}
IntelLLVM-Latest:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm-intel
    variables: {MACHINE: INTELLLVM}
PGI-21-03:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-pgi-21-03
    variables: {MACHINE: PGI}
macGNU:
    tags: ['macos']
    variables: {MACHINE: GNU}
""", yaml.Loader)


def prep_runs(test_specs, env_name, env_spec):
    jobname_compile = f'{env_name}_compile'
    pipeline_config[jobname_compile] = {
        **{'extends': '.compile_template'},
        **copy.deepcopy(env_spec),
    }

    pipeline_config[f'{env_name}_run'] = {
        **{'extends': '.simulation_template'},
        **copy.deepcopy(env_spec),
        'needs': [jobname_compile],
        'parallel': {'matrix':
                     [{'TEST_NAME': list(test_specs)}]}
    }

if __name__ == "__main__":
    try:
        specs_file = sys.argv[1]
    except IndexError:
        specs_file = "testsuite/test_branch_parameters.json"

    with open(specs_file, 'r', encoding='UTF-8') as f:
        test_specs = json.load(f)

    for env_name, env_spec in ENVIRONMENTS.items():
        prep_runs(test_specs, env_name, env_spec)

    with open('generated-config.yml', 'w', encoding='UTF-8') as f:
        f.write(yaml.dump(pipeline_config))
