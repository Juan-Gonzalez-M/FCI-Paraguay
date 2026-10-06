#!/bin/bash
# ============================================================================
# count_prose.sh — the BODY-PROSE word-count rule, saved so the figure quoted
# in the response letter is reproducible by the editor or a referee.
#
# WHY A SECOND RULE.  count_words.sh counts every token before the References
# heading in the compiled PDF, so it includes table cells, captions and table
# notes.  Reviewer 1's Comment 6 was about repetitive *prose*, and between the
# original submission and this revision the table notes grew while the prose
# shrank, so the all-inclusive figure understates the consolidation.  This
# script measures the complementary object.
#
# RULE (applied identically to every version, on the LaTeX source):
#   Take everything from \section{Introduction} to \section*{Funding}.
#   Drop every line inside a table or figure environment.
#   Drop comment lines, \begin/\end lines and displayed-equation bodies.
#   Strip LaTeX control sequences and braces, then count whitespace tokens.
#
# So: numbered-section running prose only — no tables, no captions, no table
# notes, no displayed equations, no bibliography, no front or back matter.
#
# Usage:  bash count_prose.sh [file.tex ...]
# Default: the archived V1 submission and the current submission source.
# ============================================================================
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DEFAULT_TEX=(
  "$REPO/output/archive/IREF_v1_2026-06/FCI_Paraguay_IREF_Submission.tex"
  "$REPO/output/submission/FCI_Paraguay_IREF_Submission.tex"
)

if [ "$#" -gt 0 ]; then TEXS=("$@"); else TEXS=("${DEFAULT_TEX[@]}"); fi

printf '%-58s %9s\n' "FILE" "PROSE-W"
printf '%-58s %9s\n' "----" "-------"

for f in "${TEXS[@]}"; do
  [ -f "$f" ] || { printf '%-58s %9s\n' "$(basename "$(dirname "$f")")/$(basename "$f")" "MISSING"; continue; }
  n=$(awk '
    /\\section\{Introduction\}/       { inbody = 1 }
    /\\section\*\{Funding\}/          { inbody = 0 }
    inbody != 1                       { next }
    /^[[:space:]]*%/                  { next }
    /\\begin\{(table|figure|equation)/ { skip++ ; next }
    /\\end\{(table|figure|equation)/   { if (skip > 0) skip-- ; next }
    skip > 0                          { next }
    /^[[:space:]]*\\(begin|end)\{/    { next }
    {
      line = $0
      gsub(/\\[a-zA-Z]+\*?(\[[^]]*\])?/, " ", line)   # control sequences
      gsub(/[{}$&~^_\\]/, " ", line)                  # markup punctuation
      m = split(line, w, /[[:space:]]+/)
      for (i = 1; i <= m; i++) if (w[i] ~ /[A-Za-z0-9]/) n++
    }
    END { print n + 0 }
  ' "$f")
  printf '%-58s %9s\n' "$(basename "$(dirname "$f")")/$(basename "$f")" "$n"
done
