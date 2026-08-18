#!/bin/sh
# Install CheckM2 and its UniRef100/KO reference database.
#
# CheckM2 supplies the completeness AND contamination estimates that
# completeness-aware-caller gates on. BUSCO gives completeness only, which is
# why Gotcha #4 in SKILL.md could not be acted on before this.
#
# Idempotent: re-running with the environment and database already in place
# skips both the solve and the 3.1 GB download, and still verifies the install.
#
# Usage:
#   sh scripts/install_checkm2.sh                 # env "checkm2", DB in ~/databases
#   sh scripts/install_checkm2.sh --db-path /data # custom database location
#   sh scripts/install_checkm2.sh --env-name cm2  # custom conda env name
#
# Honours an existing $CHECKM2DB: if it points at a readable .dmnd file, the
# download is skipped entirely and that database is registered.
set -eu

ENV_NAME="checkm2"
DB_PATH="${HOME}/databases"
CHECKM2_VERSION="1.1.0"
# CheckM2 1.1.0 reference database — Zenodo record 14897628.
DB_URL="https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content"
DB_TARBALL="checkm2_database.tar.gz"
DB_RELPATH="CheckM2_database/uniref100.KO.1.dmnd"

while [ $# -gt 0 ]; do
  case "$1" in
    --db-path)  DB_PATH="$2"; shift 2 ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    -h|--help)  sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

say() { printf '==> %s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- conda/mamba
if command -v mamba >/dev/null 2>&1; then
  CONDA_BIN="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BIN="conda"
else
  die "neither mamba nor conda found on PATH. Install Miniforge first:
       https://github.com/conda-forge/miniforge"
fi

# --------------------------------------------------------------- 1. environment
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  say "conda env '$ENV_NAME' already exists — skipping create"
else
  say "creating conda env '$ENV_NAME' with checkm2=$CHECKM2_VERSION"
  "$CONDA_BIN" create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda "checkm2=$CHECKM2_VERSION"
fi

INSTALLED="$(conda run -n "$ENV_NAME" checkm2 --version 2>/dev/null | tr -d '[:space:]')" \
  || die "checkm2 is not runnable inside env '$ENV_NAME'"
say "checkm2 version: $INSTALLED"

# ----------------------------------------------------------------- 2. database
# An externally supplied CHECKM2DB wins over the managed location.
if [ -n "${CHECKM2DB:-}" ] && [ -r "${CHECKM2DB:-}" ]; then
  say "using database from \$CHECKM2DB: $CHECKM2DB"
  DB_FILE="$CHECKM2DB"
else
  DB_FILE="$DB_PATH/$DB_RELPATH"
  if [ -r "$DB_FILE" ]; then
    say "database already present — skipping download ($DB_FILE)"
  else
    say "downloading reference database (~1.7 GB compressed, 3.1 GB extracted)"
    mkdir -p "$DB_PATH"
    # -C - resumes a partial download rather than restarting it.
    curl -L --retry 3 --retry-delay 5 -C - \
      -o "$DB_PATH/$DB_TARBALL" "$DB_URL"

    say "verifying archive integrity"
    gzip -t "$DB_PATH/$DB_TARBALL" \
      || die "downloaded archive is corrupt: $DB_PATH/$DB_TARBALL
       Delete it and re-run to retry the download."

    say "extracting"
    tar -xzf "$DB_PATH/$DB_TARBALL" -C "$DB_PATH"
    [ -r "$DB_FILE" ] || die "extraction did not produce $DB_FILE"

    # The tarball ships a CONTENTS.json with the expected md5; check it when
    # md5sum is available. 3 GB hash takes a few seconds and is worth it.
    if command -v md5sum >/dev/null 2>&1 && [ -r "$DB_PATH/CONTENTS.json" ]; then
      # CONTENTS.json maps the archive member to its md5; the only 32-hex
      # string in it is that digest.
      EXPECTED="$(grep -o '"[0-9a-f]\{32\}"' "$DB_PATH/CONTENTS.json" \
                  | head -1 | tr -d '"' || true)"
      if [ -n "$EXPECTED" ]; then
        say "checking md5 (this takes a moment for 3 GB)"
        ACTUAL="$(md5sum "$DB_FILE" | cut -d' ' -f1)"
        [ "$ACTUAL" = "$EXPECTED" ] \
          || die "database md5 mismatch: expected $EXPECTED, got $ACTUAL"
        say "md5 OK"
      fi
    fi

    rm -f "$DB_PATH/$DB_TARBALL"
  fi
fi

say "registering database location with CheckM2"
conda run -n "$ENV_NAME" checkm2 database --setdblocation "$DB_FILE" \
  || die "CheckM2 rejected the database at $DB_FILE"

# ------------------------------------------------------------------ 3. testrun
# This is the real gate. CheckM2 scores three known genomes and compares them
# against internal margins, so a passing testrun means the models and the
# database actually agree — not merely that the binary starts.
say "running CheckM2 testrun to verify the installation"
if conda run -n "$ENV_NAME" checkm2 testrun --threads "$(nproc 2>/dev/null || echo 4)"; then
  say "testrun passed"
else
  die "checkm2 testrun FAILED — the installation is not trustworthy.
       Do not use these scores. Check the log above for the failing stage."
fi

cat <<EOF

==> CheckM2 is ready.

    environment : $ENV_NAME
    database    : $DB_FILE

    Score genomes with:
      sh scripts/run_checkm2.sh --input <genome_dir> --output <out_dir>

    Then gate calls on the result:
      python3 completeness_aware_caller.py \\
        --checkm2-tsv <out_dir>/quality_report.tsv \\
        --genes examples/demo_genes.json --output <report_dir>
EOF
