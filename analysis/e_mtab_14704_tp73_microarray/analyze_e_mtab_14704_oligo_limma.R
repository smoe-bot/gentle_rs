#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  analysis_dir <- getwd()
}

.libPaths(c(file.path(analysis_dir, "Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(pd.clariom.d.human)
})

resource_dir <- normalizePath(file.path(
  analysis_dir,
  "../../data/publication_resources/rostock_p73_clariomd_e_mtab_14704"
))
cel_files <- list.files(resource_dir, pattern = "\\.CEL$", full.names = TRUE)
cel_files <- cel_files[order(basename(cel_files))]

sample <- data.frame(
  file = cel_files,
  sample = sub("\\.CEL$", "", basename(cel_files)),
  stringsAsFactors = FALSE
)
sample$condition <- sub("^P_SKMel29_Ad", "", sample$sample)
sample$condition <- sub("_[0-9]+$", "", sample$condition)
sample$condition <- factor(sample$condition, levels = c("GFP", "DNp73beta", "TAp73alpha"))
sample$block <- factor(sub("^.*_([0-9]+)$", "\\1", sample$sample))

message("Reading CEL files with oligo...")
raw <- read.celfiles(sample$file)
pData(raw)$condition <- sample$condition
pData(raw)$block <- sample$block
pData(raw)$sample <- sample$sample

message("RMA summarization...")
eset <- rma(raw, target = "core")
expr <- exprs(eset)
write.csv(
  data.frame(probeset_id = rownames(expr), expr, check.names = FALSE),
  file.path(analysis_dir, "oligo_rma_core_expression.csv"),
  row.names = FALSE
)

fit_limma <- function(expr, sample, paired = FALSE) {
  if (paired) {
    design <- model.matrix(~ block + condition, data = sample)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    coef_names <- c(
      DNp73beta_vs_GFP = "conditionDNp73beta",
      TAp73alpha_vs_GFP = "conditionTAp73alpha"
    )
    out <- lapply(names(coef_names), function(name) {
      tt <- topTable(fit, coef = coef_names[[name]], number = Inf, sort.by = "none")
      tt$probeset_id <- rownames(tt)
      tt$contrast <- name
      tt
    })
    ta_vs_dn <- makeContrasts(
      TAp73alpha_vs_DNp73beta = conditionTAp73alpha - conditionDNp73beta,
      levels = design
    )
    fit2 <- contrasts.fit(lmFit(expr, design), ta_vs_dn)
    fit2 <- eBayes(fit2)
    tt <- topTable(fit2, coef = "TAp73alpha_vs_DNp73beta", number = Inf, sort.by = "none")
    tt$probeset_id <- rownames(tt)
    tt$contrast <- "TAp73alpha_vs_DNp73beta"
    out[[length(out) + 1]] <- tt
    return(do.call(rbind, out))
  }

  design <- model.matrix(~ 0 + condition, data = sample)
  colnames(design) <- sub("^condition", "", colnames(design))
  contrasts <- makeContrasts(
    DNp73beta_vs_GFP = DNp73beta - GFP,
    TAp73alpha_vs_GFP = TAp73alpha - GFP,
    TAp73alpha_vs_DNp73beta = TAp73alpha - DNp73beta,
    levels = design
  )
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)
  out <- lapply(colnames(contrasts), function(name) {
    tt <- topTable(fit, coef = name, number = Inf, sort.by = "none")
    tt$probeset_id <- rownames(tt)
    tt$contrast <- name
    tt
  })
  do.call(rbind, out)
}

message("limma models...")
unpaired <- fit_limma(expr, sample, paired = FALSE)
paired <- fit_limma(expr, sample, paired = TRUE)
write.csv(unpaired, file.path(analysis_dir, "oligo_limma_unpaired_all_core.csv"), row.names = FALSE)
write.csv(paired, file.path(analysis_dir, "oligo_limma_paired_all_core.csv"), row.names = FALSE)

targets <- read.csv(file.path(analysis_dir, "target_transcript_clusters.csv"), stringsAsFactors = FALSE)
target_ids <- unique(targets$transcript_cluster_id)

subset_and_annotate <- function(tbl) {
  sub <- tbl[tbl$probeset_id %in% target_ids, ]
  merge(targets, sub, by.x = "transcript_cluster_id", by.y = "probeset_id", all.x = FALSE)
}

write.csv(
  subset_and_annotate(unpaired),
  file.path(analysis_dir, "oligo_limma_unpaired_target_transcripts.csv"),
  row.names = FALSE
)
write.csv(
  subset_and_annotate(paired),
  file.path(analysis_dir, "oligo_limma_paired_target_transcripts.csv"),
  row.names = FALSE
)

message("Done.")
