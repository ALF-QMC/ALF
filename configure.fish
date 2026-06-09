#!/usr/bin/env fish
# configure.fish - Fish shell wrapper for configure.sh
#
# Must be sourced so the exported variables take effect in your session:
#   source configure.fish MACHINE [MODE] [STAB] [options]
#   or: . configure.fish MACHINE [MODE] [STAB] [options]
#
# All arguments are forwarded verbatim to configure.sh.
# See configure.sh for the full list of MACHINE, MODE, STAB, and option values.

set -l _alf_dir (dirname (realpath (status --current-filename)))
set -l _alf_env_file (mktemp 2>/dev/null; or mktemp -t alf_env)

# Source configure.sh inside a POSIX sh subshell rooted at the ALF directory.
# Normal output (progress messages, errors) passes through to the terminal.
# On success, the resulting environment is dumped to a temp file.
sh -c 'cd "$1" && _alf_env_file="$2" && shift 2 && . ./configure.sh "$@" && env > "$_alf_env_file"' -- "$_alf_dir" "$_alf_env_file" $argv
set -l _alf_exit $status

if test $_alf_exit -ne 0
    rm -f $_alf_env_file
    return $_alf_exit
end

# Import every ALF_* variable that configure.sh exported into the fish session.
for _alf_line in (grep '^ALF_' $_alf_env_file)
    set -l _alf_kv (string split -m 1 '=' -- $_alf_line)
    if test (count $_alf_kv) -eq 2
        set -gx $_alf_kv[1] $_alf_kv[2]
    end
end

rm -f $_alf_env_file
