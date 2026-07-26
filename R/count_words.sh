#!/bin/bash
# ============================================================================
# count_words.sh — the manuscript word-count rule, saved so the figure quoted
# in the response letter is reproducible by the editor or a referee.
#
# RULE (applied identically to every version):
#   pdftotext -layout <pdf> | count every whitespace-delimited token appearing
#   BEFORE the line containing only the heading "References".
#
# This counts body prose, table cells, captions, notes and declarations, and
# excludes the bibliography.  Page counts are reported two ways because the
# two are easily confused: total PDF pages, and the page on which References
# begins (i.e. the number of pages of main text).
#
# Usage:  bash count_words.sh [pdf ...]
# Default: the archived V1 submission and the current submission PDF.
# ============================================================================
set -u

# Resolve the repository root from this script's own location, so the defaults
# work whether the script is invoked from R/, from the repository root, or by
# absolute path. Without this the relative paths silently report MISSING.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DEFAULT_PDFS=(
  "$REPO/output/archive/IREF_v1_2026-06/FCI_Paraguay_IREF_Submission.pdf"
  "$REPO/output/archive/pre_round17/FCI_Paraguay_IREF_Submission.pdf"
  "$REPO/output/submission/FCI_Paraguay_IREF_Submission.pdf"
)

if [ "$#" -gt 0 ]; then PDFS=("$@"); else PDFS=("${DEFAULT_PDFS[@]}"); fi

command -v pdftotext >/dev/null 2>&1 || {
  echo "ERROR: pdftotext (poppler) not found."; exit 1; }

printf '%-58s %9s %7s %9s\n' "FILE" "WORDS" "PAGES" "MAIN-PP"
printf '%-58s %9s %7s %9s\n' "----" "-----" "-----" "-------"

for f in "${PDFS[@]}"; do
  [ -f "$f" ] || { printf '%-58s %9s\n' "$(basename "$(dirname "$f")")/$(basename "$f")" "MISSING"; continue; }
  words=$(pdftotext -layout "$f" - 2>/dev/null | \
    awk 'BEGIN{n=0} /^[[:space:]]*References[[:space:]]*$/{exit} {n+=NF} END{print n}')
  pages=$(pdfinfo "$f" 2>/dev/null | awk '/^Pages/{print $2}')
  mainpp=$(pdftotext -layout "$f" - 2>/dev/null | \
    awk 'BEGIN{p=1} /\f/{p++} /^[[:space:]]*References[[:space:]]*$/{print p; exit}')
  printf '%-58s %9s %7s %9s\n' \
    "$(basename "$(dirname "$f")")/$(basename "$f")" "$words" "$pages" "$mainpp"
done
