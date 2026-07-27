#!/usr/bin/env Rscript
## Regenerate the bundled tissue query masks in ../data/ from the MM285 tissue
## signature. Run once; the outputs (data/<tissue>_hypo.txt and .cg) are what
## the runnable example uses, so you do NOT need R or sesame to run the example
## itself -- only to rebuild these inputs from scratch.
##
## The signature (sesameData 'MM285.tissueSignature') labels each MM285 probe
## with the tissue whose methylation it discriminates (branch) and the
## direction (type = Hypo/Hyper). A tissue's Hypo probes are the CpGs that lose
## methylation specifically in that tissue -- the marks of its regulatory
## program. We turn each Hypo set into a binary .cg over the MM285 probe
## ordering, the row space kycg tests against.
##
## Requires: R with sesame/sesameData, and `yame` + the MM285 ordering on disk
## (yame fetch MM285 -g ordering, or kycg fetch MM285:<anything>).

suppressMessages({library(sesame); library(SummarizedExperiment)})

data_dir = normalizePath(file.path(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))),
  "..", "data"), mustWork=FALSE)
if (!dir.exists(data_dir)) data_dir = "data"   # fallback: run from repo root

## the MM285 probe ordering (row space); point at your store
store = Sys.getenv("YAME_DATA_HOME",
                   file.path(Sys.getenv("HOME"), ".local/share/yame"))
ordering_gz = file.path(store, "MM285", "MM285.ordering.tsv.gz")
stopifnot(file.exists(ordering_gz))
order_ids = read.table(ordering_gz, header=TRUE, sep="\t")$Probe_ID  # 287,692

se = sesameDataGet("MM285.tissueSignature")
df = as.data.frame(rowData(se))

tissues = c("fetal_brain", "fetal_liver", "hepatocyte")
for (b in tissues) {
  q = df$Probe_ID[df$branch == b & df$type == "Hypo"]
  writeLines(q, file.path(data_dir, paste0(b, "_hypo.txt")))

  ## binary mask over the ordering: 1 where the probe is in the Hypo set
  bits = as.integer(order_ids %in% q)
  cg = file.path(data_dir, paste0(b, "_hypo.cg"))
  p = pipe(sprintf("yame pack -fb - > %s", shQuote(cg)), open="w")
  writeLines(as.character(bits), p); close(p)

  cat(sprintf("%-12s %d probes -> %s\n", b, length(q), basename(cg)))
}

## The metagene example (scripts/metagene_profile.sh) needs a beta-valued query
## instead of a binary set: one real MM285 sample's methylation, packed as a
## fraction .cg (format 4, `yame pack -fn`, NA-aware) over the same ordering.
sdf = sesameDataGet("MM285.1.SigDF")
betas = getBetas(sdf)[order_ids]                       # align to ordering
cg = file.path(data_dir, "mm285_sample.cg")
p = pipe(sprintf("yame pack -fn - > %s", shQuote(cg)), open="w")
writeLines(ifelse(is.na(betas), "NA", sprintf("%.4f", betas)), p); close(p)
cat(sprintf("mm285_sample %d probes (%d non-NA) -> %s\n",
            length(betas), sum(!is.na(betas)), basename(cg)))
