#!/usr/bin/env python3
import os
import sys
import yaml
import copy

pipeline_config = yaml.load("""
default:
  tags:
    - k8s
  artifacts:
    expire_in: 1 day

stages:
  - compile
  - simulation
  - analyze
  - summarize

summarize:
  image: git.physik.uni-wuerzburg.de:25812/z03/pdi/debian:bookworm-buildd
  stage: summarize
  when: always
  script:
    - apt-get update && apt-get install -y imagemagick
    - cd testsuite/test_vs_ed
    - ls *
    - convert */*.png -append results.png
  artifacts:
    paths:
      - testsuite/test_vs_ed
    when: always
    expire_in: 1 week

.compilation_template:
  stage: compile
  script:
    - if [ ! -O . ]; then sudo chown -R "$(id -u)" .; fi
    - export PATH="$HOME/.local/bin:$PATH"
    - pip install --no-deps pyALF
    - cd $TEST_DIR
    - ../compile.py $MACHINE
  needs:
    - pipeline: $PARENT_PIPELINE_ID
      job: prep_ED_test
  artifacts:
    paths:
      - $TEST_DIR

.simulation_template:
  stage: simulation
  script:
    - cd $TEST_DIR/$CI_NODE_INDEX
    - ../ALF.out
  artifacts:
    paths:
      - $TEST_DIR/$CI_NODE_INDEX
    exclude:
      - $TEST_DIR/$CI_NODE_INDEX/confout_*.h5
      - $TEST_DIR/ALF.out

.analysis_template:
  stage: analyze
  script:
    - export PATH="$HOME/.local/bin:$PATH"
    - pip install --no-deps pyALF
    - cd $TEST_DIR
    - ../analysis.py
  artifacts:
    paths:
      - $TEST_DIR/results.json
      - $TEST_DIR/test.png
      - $TEST_DIR/spec.yaml
    when: always
""", yaml.Loader)


def prep_runs(test_name, test_spec, env_name, env_spec):
    env_spec = copy.deepcopy(env_spec)
    test_dir = f"testsuite/test_vs_ed/{test_name}_{env_name}"
    env_spec['variables']['TEST_DIR'] = test_dir
    os.makedirs(test_dir)
    with open(os.path.join(test_dir, "spec.yaml"), 'w', encoding='UTF-8') as f:
        f.write(yaml.dump({
            'test_name': test_name,
            'env_name': env_name,
            **test_spec,
            }))

    jobname_compile = f'{test_name}_{env_name}_compile'
    pipeline_config[jobname_compile] = {
        'extends': '.compilation_template',
        **env_spec,
     }

    N_sims = len(test_spec["sim_dicts"])
    jobname_simulate=f'{test_name}_{env_name}_run'
    pipeline_config[jobname_simulate] = {
        'extends': '.simulation_template',
        **env_spec,
        'needs': [jobname_compile],
        'parallel': N_sims,
    }

    pipeline_config[f'{test_name}_{env_name}_analyze'] = {
        'extends': '.analysis_template',
        **env_spec,
        'needs': [jobname_compile, jobname_simulate],
    }


if __name__ == "__main__":
    try:
        specs_file = sys.argv[1]
    except IndexError:
        specs_file = "testsuite/test_vs_ed/test_specs.yaml"
    with open(specs_file, 'r', encoding='UTF-8') as f:
        test_specs = yaml.load(f, yaml.Loader)

    for test_name, test_spec in test_specs.items():
        for env_name, env_spec in test_spec['environments'].items():
            prep_runs(test_name, test_spec, env_name, env_spec)

    with open('generated-config.yml', 'w', encoding='UTF-8') as f:
        f.write(yaml.dump(pipeline_config))
