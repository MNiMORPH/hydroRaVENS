#!/bin/bash
# Usage: bash archive_run.sh <run-name>
#
# Copies the configuration and outputs of the current Dakota run into
# runs/<run-name>/ for version-controlled storage.
#
# Files archived:
#   dakota.in, driver.py, cannon_cfg_template.yml, run_driver.sh  -- exact config
#   evaluations.dat  (dakota.dat renamed to dodge .gitignore)      -- all evaluations
#   dakota_log.txt   (dakota.out renamed)                          -- Dakota log
#   best_fit.png     if present                                    -- diagnostic plot

set -euo pipefail

NAME="${1:?Usage: bash archive_run.sh <run-name>}"
DEST="runs/${NAME}"

if [[ -d "$DEST" ]]; then
    echo "Error: $DEST already exists. Choose a different name." >&2
    exit 1
fi

mkdir -p "$DEST"

cp dakota.in              "$DEST/"
cp driver.py              "$DEST/"
cp cannon_cfg_template.yml "$DEST/"
cp run_driver.sh          "$DEST/"
cp dakota.dat             "$DEST/evaluations.dat"
cp dakota.out             "$DEST/dakota_log.txt"
[[ -f best_fit.png ]] && cp best_fit.png "$DEST/"

N=$(( $(wc -l < "$DEST/evaluations.dat") - 1 ))
echo "Archived to $DEST  ($N evaluations)"
