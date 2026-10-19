# Independent numerical oracle using base R's predict.lm on the original models.
# Usage: Rscript scripts/export_growth_reference.R ../gRodon2/R/sysdata.rda tests/data/reference_growth_models.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Provide the source sysdata.rda and destination CSV")
original <- new.env()
load(args[1], envir = original)
features <- data.frame(
  CUBHE = c(0.8, 0.62, 0.9, 1.1, 0.5),
  ConsistencyHE = c(0.5, 0.5, 0.6, 0.55, 0.4),
  CPB = c(0.02, -0.4, -0.2, -0.1, -0.5),
  CUB = c(0.4, 0.6, 0.7, 0.65, 0.5),
  GCdiv = c(0.1, 0.05, 0.15, 0.1, 0.2),
  OGT = c(30, 25, 37, 60, 15),
  MILC = c(0.8, 0.6, 0.8, 1.2, 0.5),
  ENCprime = c(40, 40, 35, 30, 45),
  B = c(0.3, 0.2, 0.3, 0.4, 0.1),
  SCUO = c(0.4, 0.1, 0.2, 0.3, 0.1),
  MCB = c(0.2, 0.1, 0.15, 0.2, 0.07)
)
features$dCUB <- (features$CUB - features$CUBHE) / features$CUB
lambda_for <- function(name) {
  suffix <- if (grepl("_i$", name)) "_i" else ""
  if (grepl("_t150$", name)) return(original$lambda_t150)
  if (grepl("_t240$", name)) return(original$lambda_t240)
  if (grepl("_AOANOB$", name)) return(original$lambda_milc_AOANOB)
  if (grepl("_euk", name)) return(original[[paste0("lambda_milc_euk", suffix)]])
  if (grepl("_newmeta", name)) return(original[[paste0("lambda_newmeta", suffix)]])
  if (grepl("_madin", name)) return(original[[paste0("lambda_milc_madin", suffix)]])
  original$lambda_milc
}
results <- lapply(sort(ls(original, pattern = "^gRodon_model_")), function(name) {
  model <- original[[name]]
  transformed <- stats::predict.lm(model, newdata = features, interval = "confidence", level = 0.95)
  lambda <- lambda_for(name)
  hours <- (lambda * transformed + 1)^(1 / lambda)
  data.frame(model = name, case = seq_len(nrow(features)), features,
             fit = transformed[, "fit"], lower = transformed[, "lwr"], upper = transformed[, "upr"],
             hours = hours[, "fit"], hours_lower = hours[, "lwr"], hours_upper = hours[, "upr"],
             lambda = lambda, df_residual = model$df.residual,
             r_version = as.character(getRversion()), stringsAsFactors = FALSE)
})
result <- do.call(rbind, results)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(result, args[2], row.names = FALSE, na = "NaN")
cat(nrow(result), "independent R predictions written to", args[2], "\n")
