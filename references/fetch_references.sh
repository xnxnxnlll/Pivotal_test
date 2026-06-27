set -uo pipefail   # NOTE: no -e, so one failed row doesn't kill the run
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW="$HERE/raw"
MANIFEST="$HERE/manifest.tsv"
mkdir -p "$RAW"

command -v datasets >/dev/null 2>&1 || { echo "ERROR: 'datasets' not found. conda install -c conda-forge ncbi-datasets-cli" >&2; exit 1; }
command -v unzip    >/dev/null 2>&1 || { echo "ERROR: 'unzip' not found. conda install -c conda-forge unzip" >&2; exit 1; }

extract() {  # unzip a datasets bundle and append its genomic.fna to $2
  local zip="$1" dest="$2" tmp; tmp="$(mktemp -d)"
  unzip -qo "$zip" -d "$tmp" || { rm -rf "$tmp" "$zip"; return 1; }
  find "$tmp" -name '*.fna' -exec cat {} + >> "$dest"
  rm -rf "$tmp" "$zip"
}

fail=0
# strip comments/header; tr removes any CR (Windows line endings)
grep -v '^#' "$MANIFEST" | tr -d '\r' | awk -F'\t' 'NR>1 && NF>=3' \
| while IFS=$'\t' read -r taxon mode query notes; do
  # trim leading/trailing whitespace from query
  query="$(echo "$query" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
  out="$RAW/$taxon.fasta"; zip="$RAW/.${taxon}.zip"; : > "$out"
  echo ">> $taxon ($mode): $query"
  if [ "$mode" = "accession" ]; then
    # shellcheck disable=SC2086
    datasets download virus genome accession $query --include genome --filename "$zip"
  elif [ "$mode" = "taxon" ]; then
    datasets download virus genome taxon "$query" --refseq --include genome --filename "$zip"
  else
    echo "   skip: unknown mode '$mode'" >&2; continue
  fi
  if [ $? -ne 0 ] || [ ! -s "$zip" ]; then
    echo "   !! download failed for $taxon — check the accession version(s) above" >&2
    continue
  fi
  extract "$zip" "$out" && echo "   -> $(grep -c '^>' "$out") contig(s) in $(basename "$out")"
done

echo
echo "Done. If any row showed '!! download failed', fix that accession's version in"
echo "manifest.tsv and re-run. Then build the panel:  python references/build_reference.py"
