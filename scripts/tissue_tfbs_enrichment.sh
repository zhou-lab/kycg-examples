#!/usr/bin/env bash
## Tissue identity from DNA methylation — the kycg "recovers known biology"
## worked example behind the docs. A tissue's hypomethylated CpGs (from the
## MM285 tissue signature) are tested for enrichment in transcription-factor
## binding sites (ReMap). The top TFs are that tissue's lineage regulators:
## fetal brain -> neuronal factors (DLX/OTX/PAX6/LHX), hepatocyte -> hepatic
## factors (HNF4A/PPARA/RXRA/PROX1). Regenerates results/ and figures/.
##
## Tools come from PATH (conda install kycg cinderplot yame); override with the
## KYCG / CINDERPLOT / YAME env vars to point at a source build. The query
## masks in data/ are bundled; scripts/make_tissue_queries.R documents how they
## were derived from the signature (needs R+sesame, not needed to run this).
set -euo pipefail

KYCG="${KYCG:-kycg}"
CINDERPLOT="${CINDERPLOT:-cinderplot}"
YAME="${YAME:-yame}"                       # only used via `kycg fetch`

cd "$(dirname "$0")/.."                     # repo root
mkdir -p results figures

## knowledgebase: MM285 TF binding sites (ReMap), 638 factors on the MM285 row
## space. Fetched into the shared YAME store; -f skips the browser.
"$KYCG" fetch -f MM285:TFBSrm
STORE="${YAME_DATA_HOME:-$HOME/.local/share/yame}"
TFBS="$STORE/MM285/KYCG/TFBSrm.20221005.cm"

## palette
BRAIN_COL="#3b6fb6"; HEP_COL="#d1892b"     # query tissues
ENR="#c0392b"; DEP="#4a5b6b"; NS="#c8ccd0" # volcano significance

## columns of a kycg test row: 4=db(TF) 9=overlap 10=estimate(log2 OR)
##                             13=fdr 14=neglog10_fdr
for b in fetal_brain hepatocyte; do
  "$KYCG" test           -m "$TFBS" "data/${b}_hypo.cg" > "results/${b}.tsv"
  "$KYCG" test -a two.sided -m "$TFBS" "data/${b}_hypo.cg" > "results/${b}_2s.tsv"
done

## ------------------------------------------------------------------ figure 1
## Paired dot: the union of each tissue's top-8 enriched TFs, plotted for BOTH
## tissues. Each TF is enriched in one lineage and absent in the other -- the
## signatures are disjoint. (`sort | awk NR<=8` rather than `head` so pipefail
## does not trip on the closed pipe.)
top8() { awk -F'\t' 'NR>1 && $13<0.05 && $10>0 {print $14"\t"$4}' "$1" \
         | sort -k1,1gr | awk -F'\t' 'NR<=8{print $2}'; }
top8 results/fetal_brain.tsv > results/.union
top8 results/hepatocyte.tsv >> results/.union
sort -u results/.union -o results/.union

{ printf 'tf\ttissue\tneglog10_fdr\toverlap\n'
  awk -F'\t' -v t="fetal brain" 'FNR==NR{u[$1]=1;next}
       FNR>1 && ($4 in u){print $4"\t"t"\t"$14"\t"$9}' \
       results/.union results/fetal_brain.tsv
  awk -F'\t' -v t="hepatocyte" 'FNR==NR{u[$1]=1;next}
       FNR>1 && ($4 in u){print $4"\t"t"\t"$14"\t"$9}' \
       results/.union results/hepatocyte.tsv
} > results/paired_top.tsv

## level order: TF sorted by the significance of its dominant tissue; reversed
## so the top hit lands at the TOP of the coord_flip'd panel.
LV=$(awk -F'\t' 'NR>1{ if($2=="fetal brain") b[$1]=$3; else h[$1]=$3 }
       END{ for(t in b){ bb=b[t]+0; hh=h[t]+0;
            print (bb>=hh?bb:-hh)"\t"t } }' results/paired_top.tsv \
     | sort -k1,1g | awk -F'\t' '{printf "%s\"%s\"", (NR>1?",":""), $2}')

"$CINDERPLOT" "results/paired_top.tsv
  + aes(x=factor(tf, levels=c($LV)), y=neglog10_fdr, colour=tissue, size=overlap)
  + geom_point()
  + coord_flip()
  + scale_colour_manual(values=c(\"fetal brain\"=\"$BRAIN_COL\", \"hepatocyte\"=\"$HEP_COL\"))
  + labs(title=\"Tissue hypomethylation recovers lineage TFs (MM285, TFBS/ReMap)\",
         x=\"\", y=\"-log10 FDR (enrichment)\", colour=\"query tissue\")
  + theme_bw()" figures/tissue_tfbs.pdf --size 7.5x5.5
"$CINDERPLOT" "results/paired_top.tsv
  + aes(x=factor(tf, levels=c($LV)), y=neglog10_fdr, colour=tissue, size=overlap)
  + geom_point() + coord_flip()
  + scale_colour_manual(values=c(\"fetal brain\"=\"$BRAIN_COL\", \"hepatocyte\"=\"$HEP_COL\"))
  + labs(title=\"Tissue hypomethylation recovers lineage TFs (MM285, TFBS/ReMap)\",
         x=\"\", y=\"-log10 FDR (enrichment)\", colour=\"query tissue\")
  + theme_bw()" figures/tissue_tfbs.png --size 7.5x5.5

## ------------------------------------------------------------------ figure 2
## Volcano: fetal-brain hypomethylation, two-sided. Enrichment (right) names
## the top-8 neuronal TFs; the depleted cloud (left) is CpGs unmethylated in
## all tissues. |estimate|>=15 (zero-overlap saturation at log2(DBL_MIN)) is
## dropped so the finite cloud is legible.
TOP=$(awk -F'\t' 'NR>1 && $13<0.05 && $10>0 {print $14"\t"$4}' \
      results/fetal_brain_2s.tsv | sort -k1,1gr | awk -F'\t' 'NR<=8{print $2}' \
      | paste -sd'|')
awk -F'\t' -v top="$TOP" '
  BEGIN{n=split(top,a,"|"); for(i=1;i<=n;i++) T[a[i]]=1
        print "estimate\tneglog10_fdr\tsig\tlabel"}
  FNR>1{ e=$10+0; if(e<=-15||e>=15) next
         sig=($13<0.05&&e>0)?"enriched":(($13<0.05&&e<0)?"depleted":"n.s.")
         print e"\t"$14"\t"sig"\t"(($4 in T)?$4:"") }' \
  results/fetal_brain_2s.tsv > results/fetal_brain_volcano.tsv

for ext in pdf png; do
  "$CINDERPLOT" "results/fetal_brain_volcano.tsv
    + aes(estimate, neglog10_fdr, colour=sig, label=label)
    + geom_point() + geom_text()
    + scale_colour_manual(values=c(\"enriched\"=\"$ENR\", \"depleted\"=\"$DEP\", \"n.s.\"=\"$NS\"))
    + labs(title=\"Fetal-brain hypomethylation vs TF binding (MM285, two-sided)\",
           x=\"log2 odds ratio\", y=\"-log10 FDR\", colour=\"\")
    + theme_bw()" "figures/brain_volcano.$ext" --size 7x5
done

## ------------------------------------------------------------------ figure 3
## Bar: hepatocyte hypomethylation, top-10 enriched TFs (hepatic + circadian).
## A third geom; each TF on its own row, no label overlap.
awk -F'\t' 'NR>1 && $13<0.05 && $10>0 {print $14"\t"$4}' results/hepatocyte.tsv \
  | sort -k1,1gr | awk -F'\t' 'NR<=10{print $2"\t"$1}' > results/.hb
{ printf 'tf\tneglog10_fdr\n'; cat results/.hb; } > results/hep_bar.tsv
LVB=$(sort -k2,2g results/.hb | awk -F'\t' '{printf "%s\"%s\"", (NR>1?",":""), $1}')

for ext in pdf png; do
  "$CINDERPLOT" "results/hep_bar.tsv
    + aes(x=factor(tf, levels=c($LVB)), y=neglog10_fdr)
    + geom_col() + coord_flip()
    + labs(title=\"Hepatocyte hypomethylation — top enriched TFs (MM285, TFBS/ReMap)\",
           x=\"\", y=\"-log10 FDR (enrichment)\")
    + theme_bw()" "figures/hep_bar.$ext" --size 7x4.5
done

rm -f results/.union results/.hb
echo "wrote results/ and figures/{tissue_tfbs,brain_volcano,hep_bar}.{pdf,png}"
