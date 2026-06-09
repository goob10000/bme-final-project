#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SPECIES=(
  "Ecoli_K12_MG1655"
  "Yeast_S288C"
  "RaccoonDog_NyctereutesProcyonoides"
  "Mink_NeogaleVison"
  "Bat_RhinolophusAffinis"
  "Bat_RhinolophusSinicus"
)

ONLY_SPECIES=""
if [[ "${1:-}" == "--species" ]]; then
  if [[ -z "${2:-}" ]]; then
    echo "ERROR: --species requires a value" >&2
    exit 1
  fi
  ONLY_SPECIES="$2"
fi

run_one() {
  local species="$1"
  local expr_csv="$ROOT_DIR/Expression/${species}/ExpressionData.csv"
  local transcript_json="$ROOT_DIR/Genomes/${species}/gene_id_to_transcripts.json"
  local legacy_transcript_json="$ROOT_DIR/Genomes/${species}/genes.json"
  local out_dir="$ROOT_DIR/Results2/${species}"
  local out_expr="$out_dir/ExpressionData.csv"
  local out_npz="$out_dir/stuff.npz"

  echo ""
  echo "=== ${species} ==="

  if [[ ! -s "$expr_csv" ]]; then
    echo "SKIP: missing expression file: $expr_csv"
    return 0
  fi

  if [[ ! -s "$transcript_json" ]]; then
    if [[ -s "$legacy_transcript_json" ]]; then
      transcript_json="$legacy_transcript_json"
      echo "Using legacy transcript JSON: $transcript_json"
    else
      echo "SKIP: missing transcript JSON: $transcript_json"
      return 0
    fi
  fi

  mkdir -p "$out_dir"

  cp "$expr_csv" "$out_expr"
  echo "Copied expression table -> $out_expr"

  python3 "$ROOT_DIR/codon_frequency.py" --json_file "$transcript_json" -o "$out_npz"
  echo "Wrote codon matrix -> $out_npz"
}

if [[ -n "$ONLY_SPECIES" ]]; then
  run_one "$ONLY_SPECIES"
else
  for s in "${SPECIES[@]}"; do
    run_one "$s"
  done
fi

echo ""
echo "Done."
