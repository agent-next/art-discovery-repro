#!/usr/bin/env bash
# pipeline/lib/run.sh -- command runners shared by the pipeline wrappers.
#
# Nothing here parses a command string: `run` executes an argv vector and `run_sh`
# executes a LITERAL script whose variable parts arrive only as "$1".."$n".
# Data-derived values therefore never reach a shell parser. Callers set DRY_RUN.

# run CMD ARG...  -- print (dry run) or execute one command.
run() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '%s\n' "$*"
    else
        "$@"
    fi
}

# run_sh SCRIPT ARG...  -- for pipes/redirects. SCRIPT must be a literal that refers
# to its arguments as "$1".."$n". The dry run substitutes them for display only;
# execution goes through `bash -c SCRIPT bash ARG...`.
run_sh() {
    local script="$1" shown="$1" i=1 arg
    shift
    if [[ "$DRY_RUN" -eq 1 ]]; then
        for arg in "$@"; do
            shown="${shown//\"\$$i\"/\"$arg\"}"
            i=$((i + 1))
        done
        printf '%s\n' "$shown"
    else
        bash -c "$script" bash "$@"
    fi
}
