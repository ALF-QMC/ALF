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
""", yaml.Loader)


COMPILE_TEMPLATE = yaml.load("""
stage: compile
script:
  - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
  - git clone https://git.physik.uni-wuerzburg.de/ALF/pyALF.git
  - export PYTHONPATH="$PWD/pyALF:$PYTHONPATH"
  - export PATH="$PWD/pyALF/py_alf/cli:$PATH"
  - export ALF_DIR="$PWD"
  - . ./configure.sh $MACHINE noMPI HDF5 NO-INTERACTIVE
  - alf_test_branch.py
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
""", yaml.Loader)


SIMULATION_TEMPLATE = yaml.load("""
stage: simulation
script:
  - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
  - START_DIR=$PWD
  - git clone https://git.physik.uni-wuerzburg.de/ALF/pyALF.git
  - export PYTHONPATH="$PWD/pyALF:$PYTHONPATH"
  - cd ${START_DIR}/${TEST_DIR1}
  - mpiexec -n 4 ./ALF.out
  - cd ${START_DIR}/${TEST_DIR2}
  - mpiexec -n 4 ./ALF.out
  - ${START_DIR}/testsuite/compare_dirs.py
    --results ${START_DIR}/ALF_data/${TEST_NAME}.txt
    ${START_DIR}/${TEST_DIR1}
    ${START_DIR}/${TEST_DIR2}
artifacts:
  paths:
    - ${TEST_DIR1}
    - ${TEST_DIR2}
    - ALF_data/${TEST_NAME}.txt
  exclude:
    - ${TEST_DIR1}/**/confout_*.h5
    - ${TEST_DIR1}/**/ALF.out
    - ${TEST_DIR2}/**/confout_*.h5
    - ${TEST_DIR2}/**/ALF.out
  when: always
  expire_in: 1 week
""", yaml.Loader)


ANALYSIS_TEMPLATE = yaml.load("""
stage: analyze
script:
  - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
  - git clone https://git.physik.uni-wuerzburg.de/ALF/pyALF.git
  - export PYTHONPATH="$PWD/pyALF:$PYTHONPATH"
  - export PATH="$PWD/pyALF/py_alf/cli:$PATH"
  - export ALF_DIR="$PWD"
  - alf_test_branch.py \
    --sim_pars $ALF_DIR/testsuite/test_branch_parameters.json \
    --machine $MACHINE \
    --mpi \
    --branch_T $CI_COMMIT_BRANCH \
    --branch_R master \
    --devel \
    --no_prep \
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


IMAGES = yaml.load("""
Bullseye:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye
    machine: GNU
Bookworm:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm
    machine: GNU
#Intel21:
#    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-intel
#    machine: INTEL
IntelLatest:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/intel
    machine: INTEL
IntelLLVMLatest:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/intel
    machine: IntelLLVM
PGI-21-03:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-pgi-21-03
    machine: PGI
macGNU:
    machine: GNU
""", yaml.Loader)


def prep_runs(test_specs, image_name, image):
    try:
        env_spec = {'image': image['image']}
    except KeyError:
        # No image -> Assume running on MAC
        env_spec = {'tags': ['macos']}
    env_spec['variables'] = {'MACHINE': image['machine']}

    jobname_compile = f'{image_name}_compile'
    pipeline_config[jobname_compile] = {
        **copy.deepcopy(env_spec),
        **COMPILE_TEMPLATE,
    }

    analysis_needs = []
    for test_name in test_specs:
        jobname=f'{image_name}_{test_name}_run'
        env_spec['variables'] = {
            'TEST_NAME': test_name,
            'TEST_DIR1': f"ALF_data/{test_name}",
            'TEST_DIR2': f"ALF_data/{test_name}_test"}
        analysis_needs.append(jobname)
        pipeline_config[jobname] = {
            **copy.deepcopy(env_spec),
            'needs': [jobname_compile],
            **SIMULATION_TEMPLATE,
        }

    # env_spec['variables'] = {'MACHINE': image['machine']}
    # jobname_analyze = f'{image_name}_analyze'
    # pipeline_config[jobname_analyze] = {
    #     **copy.deepcopy(env_spec),
    #     'needs': analysis_needs,
    #     **ANALYSIS_TEMPLATE,
    # }


if __name__ == "__main__":
    try:
        specs_file = sys.argv[1]
    except IndexError:
        specs_file = "testsuite/test_branch_parameters.json"

    with open(specs_file, 'r', encoding='UTF-8') as f:
        test_specs = json.load(f)

    for image_name, image in IMAGES.items():
        prep_runs(test_specs, image_name, image)

    with open('generated-config.yml', 'w', encoding='UTF-8') as f:
        f.write(yaml.dump(pipeline_config))
