
rm(list = ls())
gc()

# -------------------------- Load packages & parse command line arguments --------------------------
library(data.table)
library(MASS)
library(sfsmisc)
library(optparse)  # Parse command line arguments (install if needed: install.packages("optparse"))

# Parse command line arguments
option_list <- list(
  make_option(c("-p", "--pheno"), type = "character", help = "Phenotype name", metavar = "character"),
  make_option(c("-b", "--boot_id"), type = "integer", help = "Bootstrap ID", metavar = "integer"),
  make_option(c("-d", "--data_path"), type = "character", help = "Path to data file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory", metavar = "character"),
  make_option(c("-r", "--sample_ratio"), type = "numeric", default = 0.9, help = "Sampling ratio [default: 0.9]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$pheno) || is.null(opt$boot_id) || is.null(opt$data_path) || is.null(opt$output_dir)) {
  stop("Must specify: --pheno phenotype name --boot_id bootstrap ID --data_path data path --output_dir output directory")
}

# Define parameters (simplify subsequent code)
pheno <- opt$pheno
boot_id <- opt$boot_id
data_path <- opt$data_path
output_dir <- opt$output_dir
sample_ratio <- opt$sample_ratio

# -------------------------- Data reading and preprocessing --------------------------
cat(paste0("[", Sys.time(), "] Start processing: phenotype = ", pheno, " | Bootstrap = ", boot_id, "\n"))

# Read data
dat_raw <- fread(data_path, na.strings = c("", "NA", "NaN"))
dat <- as.data.frame(dat_raw)

# Preprocessing
dat$sex <- factor(dat$sex)
dat$age_at_ados <- as.numeric(dat$age_at_ados)
# New: handle ados_module type (adjust according to actual data type; if numeric, change to as.numeric)
dat$ados_module <- factor(dat$ados_module)  
colnames(dat) <- gsub("-", "_", colnames(dat))
genes <- colnames(dat)[2:15534]  # Gene column range (adjust as needed)

# Filter missing phenotype + missing covariates (critical: avoid model errors)
# Construct dynamic filtering condition: for adi_r phenotypes additionally require non‑missing ados_module
filter_cond <- !is.na(dat[[pheno]]) & !is.na(dat$sex)
if (grepl("^adi_r", pheno)) {  # Check if phenotype starts with adi_r
  filter_cond <- filter_cond & !is.na(dat$age_at_ados) & !is.na(dat$ados_module)
} else {
  filter_cond <- filter_cond & !is.na(dat$age_at_ados)
}
dat_pheno <- dat[filter_cond, ]
cat(paste0("[", Sys.time(), "] Sample size after filtering: ", nrow(dat_pheno), "\n"))

# -------------------------- Core functions --------------------------
# Generate bootstrap sample
#load(paste0("../bin/bootstrap_2023_samples_list_",pheno,".Rdat"))
#load("bootstrap_samples_list.Rdat")
load("bootstrap_1711_samples_list_5_pheno.Rdat")
dat_boot_ids <- sampled_id_list[[boot_id]]
#dat_boot_ids <- random_samples[[boot_id]]
dat_boot <- dat_pheno[match(dat_boot_ids, dat_pheno$V1), ]
cat(paste0("[", Sys.time(), "] Bootstrap sample size: ", nrow(dat_boot), "\n"))

# 2. RLM analysis (with tryCatch error handling)
run_rlm <- function(gene, pheno, data) {
  # Filter out invalid gene values
  idx_valid <- !is.infinite(data[[gene]]) & !is.na(data[[gene]])
  data_valid <- data[idx_valid, ]
  
  # Core modification: dynamically construct formula based on phenotype prefix
  if (grepl("^adi_r", pheno)) {
    # adi_r phenotype: covariates = interaction of ados_module and age_at_ados + sex
    formula_str <- as.formula(paste(pheno, "~", gene, "+ ados_module:age_at_ados + sex"))
  } else {
    # Non‑adi_r phenotype: keep original covariates (age_at_ados + sex)
    formula_str <- as.formula(paste(pheno, "~", gene, "+ age_at_ados + sex"))
  }
  
  rlm_model <- rlm(formula_str, data = data_valid)
  
  # Extract coefficient results
  coef_summary <- summary(rlm_model)$coefficients
  
  # Wald test (f.robftest)
  wald_test <- f.robftest(rlm_model, var = gene)
  
  # Compile results
  res <- data.frame(
    gene = gene,
    phenotype = pheno,
    value = coef_summary[2, "Value"],
    SE = coef_summary[2, "Std. Error"],
    t_value = coef_summary[2, "t value"],
    p_value = as.numeric(wald_test$p.value),    # Wald test p‑value
    bootstrap_id = NA                           # Will be filled with bootstrap ID later
  )
  return(res)
}

# -------------------------- Execute analysis --------------------------
# Batch run RLM
rlm_res <- lapply(genes, function(g) run_rlm(g, pheno, dat_boot))
rlm_res_df <- do.call(rbind, rlm_res)

# Fill bootstrap ID
rlm_res_df$bootstrap_id <- boot_id

# -------------------------- Save results --------------------------
# Output file name: phenotype_BootstrapID.Rdat
output_file <- file.path(output_dir, sprintf("RLM_%s_boot_%d.Rdat", pheno, boot_id))
save(rlm_res_df, file = output_file)

cat(paste0("[", Sys.time(), "] Done! Results saved to: ", output_file, "\n"))
cat(paste0("Number of valid results: ", sum(!is.na(rlm_res_df$p_value)), "/", nrow(rlm_res_df), "\n"))

# Clean memory
rm(list = ls())
gc()