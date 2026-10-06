#!/bin/bash
# ============================================================================
# check_overlaps.sh — production QA: find text that physically overlaps in the
# compiled PDFs.
#
# WHY.  Pandoc sets pipe-table column widths in proportion to the dash counts
# in the separator row, with no regard for what the cells contain.  A column
# holding wide or unbreakable content (inline math such as $[-4.36,\,-0.66]$,
# or a \texttt file path) therefore overflows silently into its neighbour.  The
# result is legible in the source and broken in the PDF.
#
# Do NOT use `pdftotext` text output to look for this.  It gives false results
# in BOTH directions on this document: it splits `(0.008)` and joins
# "are preselected" where the PDF is correct (bold and \emph font switches are
# emitted as separate flow blocks), and it inserts spaces where the PDF really
# does overlap.  Glyph bounding boxes are the only reliable test.
#
# RULE: flag any two horizontally adjacent words on the same line whose
# bounding boxes overlap by more than 0.5 pt.
#
# EXPECTED NON-ZERO RESULTS.  Combining math accents (i + U+0302), and sub- or
# superscripts, legitimately share horizontal space with their base glyph.
# These are reported separately and are not defects.
#
# Usage:  bash check_overlaps.sh [file.pdf ...]
# Default: the four files of the submission package.
# ============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DEFAULT_PDFS=(
  "$REPO/output/submission/FCI_Paraguay_IREF_Submission.pdf"
  "$REPO/output/pdf/FCI_Paraguay_Online_Appendix.pdf"
  "$REPO/private/Response_to_Reviewers.pdf"
  "$REPO/output/submission/Cover_Letter_IREF.pdf"
)

if [ "$#" -gt 0 ]; then PDFS=("$@"); else PDFS=("${DEFAULT_PDFS[@]}"); fi

command -v pdftotext >/dev/null 2>&1 || {
  echo "ERROR: pdftotext (poppler) not found."; exit 1; }

TMP="$(mktemp -t overlapxml)"
trap 'rm -f "$TMP"' EXIT

for f in "${PDFS[@]}"; do
  [ -f "$f" ] || { printf '%-46s MISSING\n' "$(basename "$f")"; continue; }
  pdftotext -bbox-layout "$f" "$TMP" 2>/dev/null
  printf '%-46s ' "$(basename "$f")"
  python3 - "$TMP" <<'PY'
import re, html, sys
s = open(sys.argv[1], encoding='utf8').read()
hits = []
for pi, pg in enumerate(re.split(r'<page ', s)[1:], start=1):
    for line in re.findall(r'<line[^>]*>(.*?)</line>', pg, re.S):
        ws = re.findall(
            r'<word xMin="([\d.]+)" yMin="[\d.]+" xMax="([\d.]+)" yMax="[\d.]+">(.*?)</word>',
            line, re.S)
        for a, b in zip(ws, ws[1:]):
            gap = float(b[0]) - float(a[1])
            if gap < -0.5:
                hits.append((pi, round(gap, 2), html.unescape(a[2]), html.unescape(b[2])))

def is_accent(t):
    # combining marks, and the sub/superscript tokens this document produces
    return (any(0x300 <= ord(c) <= 0x36F for c in t)
            or t in ('\U0001d43c\U0001d447', '∗', 'placebo', '\U0001d441\U0001d44a'))

genuine = [h for h in hits if not is_accent(h[3])]
print(f"{len(hits):>3} overlaps | {len(hits)-len(genuine):>3} accent/superscript | "
      f"{len(genuine):>3} GENUINE")
for h in genuine:
    print(f"      p{h[0]:>3}  {h[1]:>7}pt   '{h[2][:34]}' || '{h[3][:34]}'")
PY
done
