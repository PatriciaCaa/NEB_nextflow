#!/usr/bin/env bash

# Record software versions and source-code checksums for reproducibility.
# This script uses CSD3-specific environment paths and requires adaptation before running on another computing system.
set -euo pipefail

OUT="${1:?Provide the provenance output directory}"
mkdir -p "$OUT/software"

BASE=/rds-d6/user/vpc25/hpc-work
MAMBA="$BASE/tools/micromamba/bin/micromamba"

# Record the installed packages in each RNA environment.
for ENV in rna rna_report alignment neb-pychopper; do
    MAMBA_ROOT_PREFIX="$BASE" "$MAMBA" list \
        -p "$BASE/envs/$ENV" --json \
        > "$OUT/software/${ENV}_packages.json"
done

# Record the patched Pychopper source checksum.
sha256sum \
    "$BASE/envs/neb-pychopper/lib/python3.13/site-packages/pychopper/scripts/pychopper.py" \
    > "$OUT/software/pychopper_source.sha256"

# Record checksums of pipeline source code and configurations.
find test_rna.nf nextflow.config conf modules/rna \
     modules/dna/provenance subworkflows scripts envs \
    -type f \( -name '*.nf' -o -name '*.config' \
        -o -name '*.R' -o -name '*.py' \
        -o -name '*.yml' -o -name '*.sh' \) \
    -print0 | sort -z | xargs -0 -r sha256sum \
    > "$OUT/code_checksums.sha256"

# Record Nextflow version.
nextflow -version > "$OUT/software/nextflow_version.txt" 2>&1

echo "Software and code provenance captured: $OUT"
