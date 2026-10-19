# Run installed original gRodon. All values and paths arrive as JSON data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide exactly one JSON parameter-file path")
required <- c("gRodon", "Biostrings", "jsonlite")
for (package in required) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(paste("Install the original R package/dependency before comparison:", package))
  }
}
parameters <- jsonlite::fromJSON(args[[1]], simplifyVector = TRUE)
genes <- Biostrings::readDNAStringSet(parameters$cds_path)
gene_ids <- sub("[[:space:]].*$", "", names(genes))
he_ids <- trimws(readLines(parameters$ids_path, warn = FALSE, encoding = "UTF-8"))
he_ids <- he_ids[nzchar(he_ids)]
if (length(he_ids) == 0L || anyDuplicated(he_ids)) stop("HE IDs must be nonempty and unique")
if (anyDuplicated(gene_ids)) stop("CDS FASTA identifiers must be unique")
if (!all(he_ids %in% gene_ids)) stop("Some HE IDs are absent from the CDS FASTA")
he <- gene_ids %in% he_ids
temperature <- if (is.null(parameters$temperature)) "none" else parameters$temperature
warnings_seen <- character()
start <- proc.time()[["elapsed"]]
incidental_output <- capture.output({
  result <- withCallingHandlers(
    gRodon::predictGrowth(
      genes, he, mode = parameters$mode, training_set = parameters$training_set,
      temperature = temperature, genetic_code = parameters$genetic_code,
      fragments = FALSE
    ),
    warning = function(warning) {
      warnings_seen <<- c(warnings_seen, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
})
elapsed <- proc.time()[["elapsed"]] - start
if (length(incidental_output)) cat(paste(incidental_output, collapse = "\n"), "\n", file = stderr())
statistics <- lapply(result, function(value) {
  if (!is.numeric(value) || length(value) != 1L) stop("gRodon returned a non-scalar statistic")
  unname(as.numeric(value))
})
versions <- list(R = R.version.string)
for (package in c("gRodon", "coRdon", "Biostrings", "matrixStats", "jsonlite")) {
  versions[[package]] <- as.character(utils::packageVersion(package))
}
out <- list(
  implementation = "original R package", versions = versions,
  settings = list(mode = parameters$mode, training_set = parameters$training_set,
                  temperature = temperature, fragments = FALSE,
                  genetic_code = parameters$genetic_code, confidence = 0.95),
  input_genes = length(genes), highly_expressed_before_filtering = sum(he),
  computation_seconds = elapsed, statistics = statistics, warnings = warnings_seen,
  session_info = capture.output(utils::sessionInfo())
)
cat(jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE, digits = NA,
                    na = "null", null = "null"), "\n")
