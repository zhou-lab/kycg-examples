#!/usr/bin/env bash
## ChromHMM enrichment of a single-cell methylome — the kycg + cinderplot
## worked example behind the docs. Regenerates results/ and figures/ from
## data/. Tools come from PATH (conda install kycg cinderplot); override with
## the KYCG / CINDERPLOT env vars to point at a source build.
set -euo pipefail

KYCG="${KYCG:-kycg}"
CINDERPLOT="${CINDERPLOT:-cinderplot}"

## run from the repo root regardless of where this is invoked from
cd "$(dirname "$0")/.."
mkdir -p results figures

## 1. enrichment: one query cell against the 18 ChromHMM states -----------
##    (query .cg and knowledgebase .cm must index the same row space)
"$KYCG" test -m data/chromhmm.cm data/onecell.cg > results/res.tsv

## 2. plot prep -----------------------------------------------------------
## A zero/inf odds ratio saturates to DBL_MIN/DBL_MAX, so estimate=log2(OR)
## lands at +-~1022 and squashes the axes; drop |estimate|>=1000 as
## KYCG_plotVolcano does. cinderplot has no inline ifelse, so precompute a
## significance column (FDR<0.05  =>  neglog10_fdr>1.30103). Columns: 10=estimate
## 14=neglog10_fdr.
awk -F'\t' 'NR==1 { print $0"\tsig"; next }
            ($10 < 1000 && $10 > -1000) {
              print $0"\t" ($14 >= 1.30103 ? "Significant" : "NotSig") }' \
    results/res.tsv > results/res_aug.tsv

## waterfall needs a rank index: rows sorted by estimate (col 10), numbered
{ head -1 results/res_aug.tsv | awk -F'\t' '{ print "index\t"$0 }'
  tail -n +2 results/res_aug.tsv | sort -t$'\t' -k10,10g |
    awk -F'\t' '{ print NR"\t"$0 }'; } > results/res_wf.tsv

## 3. figures (PDF for print, PNG for the web docs) -----------------------
## cinderplot does NOT tilde-expand inside the quoted DSL; relative paths are
## fine because we cd'd to the repo root above.
for ext in pdf png; do
  ## volcano: log2(OR) vs -log10 FDR, coloured by significance
  "$CINDERPLOT" "results/res_aug.tsv
      + aes(estimate, neglog10_fdr, colour=sig)
      + geom_point()
      + scale_colour_manual(values=c(Significant=\"#c0392b\", NotSig=\"#7f8c8d\"))
      + labs(title=\"ChromHMM enrichment (onecell)\",
             x=\"log2(OR)\", y=\"-log10 FDR\")" \
      -o "figures/volcano.$ext"

  ## dot: state vs -log10 FDR, size=overlap, colour=estimate, horizontal
  "$CINDERPLOT" "results/res_aug.tsv
      + aes(db, neglog10_fdr, size=overlap, colour=estimate)
      + geom_point() + coord_flip()
      + labs(title=\"ChromHMM enrichment (onecell)\", y=\"-log10 FDR\", x=\"\")" \
      -o "figures/dot.$ext"

  ## bar: state vs -log10 FDR, horizontal
  "$CINDERPLOT" "results/res_aug.tsv
      + aes(db, neglog10_fdr) + geom_col() + coord_flip()
      + labs(title=\"ChromHMM enrichment (onecell)\", y=\"-log10 FDR\", x=\"\")" \
      -o "figures/bar.$ext"

  ## waterfall: enrichment ranked, point size = -log10 FDR
  "$CINDERPLOT" "results/res_wf.tsv
      + aes(index, estimate, size=neglog10_fdr)
      + geom_point() + geom_hline(yintercept=0)
      + labs(title=\"ChromHMM waterfall (onecell)\", x=\"rank\", y=\"log2(OR)\")" \
      -o "figures/waterfall.$ext"
done

echo "wrote results/res.tsv and figures/{volcano,dot,bar,waterfall}.{pdf,png}"
