#!/bin/bash

set -veo pipefail

eval ${_ARRIBA_STANDALONE_INIT}

bam_input_base=$(basename ${_ARRIBA_STANDALONE_BAM} .out.bam)
bam_input_temp=${TMPDIR}/${bam_input_base}.out.bam
bam_output_temp=${TMPDIR}/${bam_input_base}.sortedByCoord.out.bam

cp ${_ARRIBA_STANDALONE_BAM} ${bam_input_temp}

samtools sort \
    -@ ${_ARRIBA_STANDALONE_SORT_THREADS} \
    -o "${bam_output_temp}" \
    "${bam_input_temp}"

samtools index \
    -@ ${_ARRIBA_STANDALONE_INDEX_THREADS} \
    "${bam_output_temp}"

mkdir -p $(dirname ${_ARRIBA_STANDALONE_OUTPUT})
draw_fusions.R \
    --fusions="${_ARRIBA_STANDALONE_FUSIONS}" \
    --alignments="${bam_output_temp}" \
    --output="${_ARRIBA_STANDALONE_OUTPUT}" \
    --annotation="${_ARRIBA_STANDALONE_ANNOTATION}" \
    --cytobands="${_ARRIBA_STANDALONE_CYTOBANDS}" \
    --proteinDomains="${_ARRIBA_STANDALONE_PROTEIN_DOMAINS}"
