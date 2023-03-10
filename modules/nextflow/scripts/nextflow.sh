#!/bin/bash

set -e -o pipefail

_clean () {
    echo "Killing nextflow..."
    kill -TERM $_nf_pid
    wait $_nf_pid
    code=$?
    exit $code
}

module load $_JAVA_MODULE
module load $_NEXTFLOW_MODULE
NXF_HOME="${TMPDIR}/.nextflow" nextflow $@ & _nf_pid=$!

trap _clean EXIT

wait $_nf_pid
exit $?
