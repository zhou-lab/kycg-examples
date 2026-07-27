#!/usr/bin/env bash
## Metagene methylation profile — kycg reads out the classic gene methylation
## architecture. A real MM285 sample's betas are averaged over the 30 metagene
## bins of the MetagenePC knowledgebase (10 upstream, the gene body in deciles,
## 10 downstream). The profile is textbook: a hypomethylation dip at the
## promoter/TSS, rising methylation across the gene body toward the 3' end, and
## a return to the intergenic baseline downstream.
##
## Tools from PATH (conda install kycg cinderplot yame); override via the
## KYCG / CINDERPLOT env vars. The beta query mask (data/mm285_sample.cg) is
## bundled; scripts/make_tissue_queries.R records how it was built (needs R +
## sesame, not needed to run this).
set -euo pipefail

KYCG="${KYCG:-kycg}"
CINDERPLOT="${CINDERPLOT:-cinderplot}"

cd "$(dirname "$0")/.."                     # repo root
mkdir -p results figures

## knowledgebase: MM285 metagene coordinates, fetched into the shared store.
"$KYCG" fetch -f MM285:MetagenePC
STORE="${YAME_DATA_HOME:-$HOME/.local/share/yame}"
MG="$STORE/MM285/KYCG/MetagenePC.20220911.cm"

## kycg test averages the sample's beta over each metagene bin (the `beta`
## column). The bin name is "<index>;<span>"; the leading index orders the
## x-axis: -10..-1 upstream, 0..9 gene body deciles, 10..19 downstream.
"$KYCG" test -m "$MG" data/mm285_sample.cg > results/metagene_raw.tsv

## columns: 4 = db (bin name), 20 = beta (mean methylation in the bin)
{ printf 'pos\tbeta\n'
  awk -F'\t' 'NR>1{split($4,a,";"); print a[1]"\t"$20}' results/metagene_raw.tsv \
    | sort -k1,1n; } > results/metagene.tsv

## line profile; the two rules mark the TSS (between -1 and 0) and the TES
## (between 9 and 10), bracketing the gene body.
for ext in pdf png; do
  "$CINDERPLOT" "results/metagene.tsv
    + aes(pos, beta)
    + geom_vline(xintercept=-0.5) + geom_vline(xintercept=9.5)
    + geom_line() + geom_point()
    + labs(title=\"Metagene methylation profile (MM285 sample)\",
           x=\"metagene position  ( -10..-1 upstream | 0..9 gene body | 10..19 downstream )\",
           y=\"mean beta\")
    + theme_bw()" "figures/metagene.$ext" --size 7.5x4.2
done

echo "wrote results/metagene.tsv and figures/metagene.{pdf,png}"
