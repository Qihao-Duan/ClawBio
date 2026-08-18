#!/bin/sh
# Score genomes with CheckM2 and emit a quality_report.tsv that
# completeness_aware_caller.py can consume via --checkm2-tsv.
#
# Usage:
#   sh scripts/run_checkm2.sh --input <genome_dir> --output <out_dir>
#   sh scripts/run_checkm2.sh --input bins/ --output out/ --threads 8 --lowmem
#
# Options:
#   --input     directory of genome FASTAs (required)
#   --output    output directory (required)
#   --threads   CPUs to use (default: all)
#   --extension FASTA extension to pick up (default: .fna)
#   --env-name  conda env holding CheckM2 (default: checkm2)
#   --lowmem    shrink the DIAMOND block size; slower but far less RAM
#
# Run scripts/install_checkm2.sh first. The reference database is NOT fetched
# here: a 3 GB download in the middle of an analysis is slow, surprising, and
# bad for reproducibility, so a missing database is an error with a pointer.
set -eu

ENV_NAME="checkm2"
INPUT=""
OUTPUT=""
THREADS="$(nproc 2>/dev/null || echo 4)"
EXTENSION=".fna"
LOWMEM=""

while [ $# -gt 0 ]; do
  case "$1" in
    --input)     INPUT="$2"; shift 2 ;;
    --output)    OUTPUT="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --extension) EXTENSION="$2"; shift 2 ;;
    --env-name)  ENV_NAME="$2"; shift 2 ;;
    --lowmem)    LOWMEM="--lowmem"; shift ;;
    -h|--help)   sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

[ -n "$INPUT" ]  || die "--input is required"
[ -n "$OUTPUT" ] || die "--output is required"
[ -d "$INPUT" ]  || die "input directory does not exist: $INPUT"

command -v conda >/dev/null 2>&1 \
  || die "conda not found on PATH. Run scripts/install_checkm2.sh first."

conda run -n "$ENV_NAME" checkm2 --version >/dev/null 2>&1 \
  || die "CheckM2 is not available in conda env '$ENV_NAME'.
       Run: sh scripts/install_checkm2.sh"

# Fail early and specifically on a missing database, rather than letting
# CheckM2 get all the way through Prodigal before discovering it.
if [ -n "${CHECKM2DB:-}" ]; then
  [ -r "$CHECKM2DB" ] \
    || die "\$CHECKM2DB is set but not readable: $CHECKM2DB"
elif ! conda run -n "$ENV_NAME" checkm2 database --current >/dev/null 2>&1; then
  die "CheckM2 reference database not found.
       Run: sh scripts/install_checkm2.sh
       or:  export CHECKM2DB=/path/to/CheckM2_database/uniref100.KO.1.dmnd"
fi

COUNT="$(find "$INPUT" -maxdepth 1 -name "*$EXTENSION" | wc -l | tr -d ' ')"
[ "$COUNT" -gt 0 ] \
  || die "no *$EXTENSION files in $INPUT (use --extension to change this)"

printf '==> scoring %s genome(s) from %s with %s threads\n' \
  "$COUNT" "$INPUT" "$THREADS"

# Intermediates are kept so --resume can reuse Prodigal and DIAMOND work if
# the run is interrupted; DIAMOND is by far the expensive stage.
conda run --live-stream -n "$ENV_NAME" checkm2 predict \
  --input "$INPUT" \
  --output-directory "$OUTPUT" \
  --threads "$THREADS" \
  --extension "$EXTENSION" \
  --force \
  $LOWMEM

REPORT="$OUTPUT/quality_report.tsv"
[ -r "$REPORT" ] || die "CheckM2 finished but produced no $REPORT"

printf '\n==> quality report: %s\n' "$REPORT"
printf '    gate calls on it with:\n'
printf '      python3 completeness_aware_caller.py --checkm2-tsv %s \\\n' "$REPORT"
printf '        --genome <name> --genes <genes.json> --output <report_dir>\n'
