#!/usr/bin/env Rscript

# =============================================================
# Exome Data Analysis Workflow (R / RStudio Script)
# =============================================================
# This script provides a practical template for exome sequencing
# variant analysis in R. It includes:
#   1) VCF import
#   2) Basic sample/variant QC summaries
#   3) Filtering by quality and consequence
#   4) Functional annotation via Ensembl VEP REST API (optional)
#   5) Case-control burden test at the gene level
#   6) Export of filtered variants and burden results
#
# Adapt paths and thresholds below for your own project.
# =============================================================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(httr2)
  library(jsonlite)
})

# -------------------------
# User configuration
# -------------------------
vcf_file <- "data/exome_cohort.vcf.gz"
output_dir <- "results"
sample_metadata_file <- "data/sample_metadata.csv" # must include: sample_id, phenotype (case/control)

# Filtering thresholds
min_dp <- 10
min_gq <- 20
max_missing_rate <- 0.1
rare_af_cutoff <- 0.01

# Optional: if your VCF has consequence and AF annotations in INFO
info_af_field <- "AF"
info_consequence_field <- "CSQ"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------
# Utility functions
# -------------------------
message_step <- function(msg) {
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

safe_read_metadata <- function(path) {
  if (!file.exists(path)) {
    stop("Sample metadata file not found: ", path)
  }

  meta <- fread(path) %>% as.data.frame()
  required_cols <- c("sample_id", "phenotype")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  meta$phenotype <- tolower(meta$phenotype)
  if (!all(meta$phenotype %in% c("case", "control"))) {
    stop("phenotype column must contain only 'case' or 'control'")
  }

  meta
}

calc_missing_rate <- function(gt_matrix) {
  # Missing genotype coded as NA or ./.
  is_missing <- is.na(gt_matrix) | gt_matrix == "./." | gt_matrix == ".|."
  rowMeans(is_missing)
}

split_gene_from_csq <- function(csq_string) {
  # Very lightweight parser for VEP CSQ annotations where gene symbol is often in field 4 or 5
  # depending on your VCF header. Adjust as needed for your VEP format.
  ifelse(
    is.na(csq_string),
    NA_character_,
    vapply(strsplit(csq_string, ","), function(x) {
      first_ann <- strsplit(x[[1]], "\\|")[[1]]
      idx <- which(first_ann != "")[1]
      if (length(idx) == 0) return(NA_character_)
      first_ann[min(5, length(first_ann))]
    }, character(1))
  )
}

# Optional: annotate one variant with Ensembl VEP REST endpoint
vep_annotate_variant <- function(chr, pos, ref, alt, species = "human") {
  hgvs <- sprintf("%s:g.%s%s>%s", chr, pos, ref, alt)
  url <- sprintf("https://rest.ensembl.org/vep/%s/hgvs/%s", species, URLencode(hgvs, reserved = TRUE))

  resp <- request(url) %>%
    req_headers("Content-Type" = "application/json", "Accept" = "application/json") %>%
    req_perform()

  if (resp_status(resp) >= 300) return(NULL)
  resp_body_json(resp, simplifyVector = TRUE)
}

# -------------------------
# 1) Import data
# -------------------------
message_step("Reading metadata")
metadata <- safe_read_metadata(sample_metadata_file)

message_step("Loading VCF")
if (!file.exists(vcf_file)) {
  stop("VCF file not found: ", vcf_file)
}
vcf <- readVcf(vcf_file, genome = "hg38")

# Ensure sample order matches metadata
vcf_samples <- samples(header(vcf))
metadata <- metadata %>% filter(sample_id %in% vcf_samples)
metadata <- metadata[match(vcf_samples, metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) {
  stop("Some VCF samples are missing in metadata file")
}

# -------------------------
# 2) QC summaries
# -------------------------
message_step("Computing variant/sample QC summaries")
gt <- geno(vcf)$GT
dp <- geno(vcf)$DP
gq <- geno(vcf)$GQ

missing_rate <- calc_missing_rate(gt)
mean_dp <- rowMeans(dp, na.rm = TRUE)
mean_gq <- rowMeans(gq, na.rm = TRUE)

qc_table <- data.frame(
  chr = as.character(seqnames(rowRanges(vcf))),
  pos = start(rowRanges(vcf)),
  ref = as.character(ref(vcf)),
  alt = sapply(alt(vcf), function(a) as.character(a[[1]])),
  missing_rate = missing_rate,
  mean_dp = mean_dp,
  mean_gq = mean_gq,
  stringsAsFactors = FALSE
)

fwrite(qc_table, file.path(output_dir, "variant_qc_metrics.tsv"), sep = "\t")

# -------------------------
# 3) Variant filtering
# -------------------------
message_step("Applying filters")
pass_missing <- missing_rate <= max_missing_rate
pass_dp <- mean_dp >= min_dp
pass_gq <- mean_gq >= min_gq

af_vec <- info(vcf)[[info_af_field]]
if (is.null(af_vec)) {
  warning("INFO/", info_af_field, " not found. Rare AF filtering skipped.")
  pass_af <- rep(TRUE, nrow(vcf))
} else {
  if (is.matrix(af_vec)) af_vec <- af_vec[, 1]
  pass_af <- af_vec <= rare_af_cutoff
  pass_af[is.na(pass_af)] <- TRUE
}

csq_vec <- info(vcf)[[info_consequence_field]]
if (is.null(csq_vec)) {
  warning("INFO/", info_consequence_field, " not found. Consequence filtering skipped.")
  pass_consequence <- rep(TRUE, nrow(vcf))
} else {
  keep_terms <- c("missense_variant", "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant")
  pass_consequence <- vapply(csq_vec, function(x) {
    any(grepl(paste(keep_terms, collapse = "|"), x, ignore.case = TRUE))
  }, logical(1))
}

keep_idx <- pass_missing & pass_dp & pass_gq & pass_af & pass_consequence
vcf_filt <- vcf[keep_idx, ]

message_step(sprintf("Retained %d / %d variants", nrow(vcf_filt), nrow(vcf)))

# -------------------------
# 4) Build variant result table
# -------------------------
message_step("Building filtered variant table")
filtered_tbl <- data.frame(
  chr = as.character(seqnames(rowRanges(vcf_filt))),
  pos = start(rowRanges(vcf_filt)),
  ref = as.character(ref(vcf_filt)),
  alt = sapply(alt(vcf_filt), function(a) as.character(a[[1]])),
  af = {
    af_f <- info(vcf_filt)[[info_af_field]]
    if (is.null(af_f)) NA_real_ else if (is.matrix(af_f)) af_f[, 1] else af_f
  },
  consequence = {
    csq_f <- info(vcf_filt)[[info_consequence_field]]
    if (is.null(csq_f)) NA_character_ else csq_f
  },
  stringsAsFactors = FALSE
)

filtered_tbl$gene <- split_gene_from_csq(filtered_tbl$consequence)
fwrite(filtered_tbl, file.path(output_dir, "filtered_variants.tsv"), sep = "\t")

# -------------------------
# 5) Gene-level burden test (case/control)
# -------------------------
message_step("Running simple gene-level burden test")

# Convert GT to carrier matrix: 1 if genotype is not 0/0 and not missing
carrier_matrix <- geno(vcf_filt)$GT
carrier_matrix <- apply(carrier_matrix, c(1, 2), function(x) {
  if (is.na(x) || x %in% c("./.", ".|.")) return(NA_integer_)
  if (x %in% c("0/0", "0|0")) return(0L)
  return(1L)
})

if (!is.matrix(carrier_matrix)) {
  carrier_matrix <- matrix(carrier_matrix, nrow = nrow(vcf_filt), byrow = FALSE)
}

gene_vec <- filtered_tbl$gene
gene_vec[is.na(gene_vec) | gene_vec == ""] <- "UNKNOWN_GENE"

is_case <- metadata$phenotype == "case"

burden_results <- lapply(unique(gene_vec), function(g) {
  idx <- which(gene_vec == g)
  gene_carrier <- colSums(carrier_matrix[idx, , drop = FALSE], na.rm = TRUE) > 0

  case_carrier <- sum(gene_carrier[is_case], na.rm = TRUE)
  case_noncarrier <- sum(is_case, na.rm = TRUE) - case_carrier
  ctrl_carrier <- sum(gene_carrier[!is_case], na.rm = TRUE)
  ctrl_noncarrier <- sum(!is_case, na.rm = TRUE) - ctrl_carrier

  cont_tbl <- matrix(c(case_carrier, case_noncarrier, ctrl_carrier, ctrl_noncarrier), nrow = 2)
  ft <- fisher.test(cont_tbl)

  data.frame(
    gene = g,
    case_carrier = case_carrier,
    case_noncarrier = case_noncarrier,
    control_carrier = ctrl_carrier,
    control_noncarrier = ctrl_noncarrier,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
})

burden_df <- bind_rows(burden_results) %>% arrange(p_value)
burden_df$fdr <- p.adjust(burden_df$p_value, method = "BH")

fwrite(burden_df, file.path(output_dir, "gene_burden_results.tsv"), sep = "\t")

# -------------------------
# 6) Export filtered VCF
# -------------------------
message_step("Writing filtered VCF")
out_vcf <- file.path(output_dir, "filtered_exome_variants.vcf.gz")
writeVcf(vcf_filt, out_vcf, index = TRUE)

message_step("Analysis complete.")
message("Output files:")
message(" - ", file.path(output_dir, "variant_qc_metrics.tsv"))
message(" - ", file.path(output_dir, "filtered_variants.tsv"))
message(" - ", file.path(output_dir, "gene_burden_results.tsv"))
message(" - ", out_vcf)
