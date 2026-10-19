# Generate independent numerical fixtures with the original gRodon R functions.
# Usage: Rscript scripts/generate_codon_reference.R /path/to/gRodon2 [output.json]
# Required R packages: coRdon, Biostrings, dplyr, magrittr, matrixStats, jsonlite.
# All sequences below are artificial and generated solely for numerical tests.

suppressPackageStartupMessages({
  library(coRdon)
  library(Biostrings)
  library(magrittr)
  library(matrixStats)
  library(dplyr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
upstream <- if (length(args) >= 1) args[[1]] else "../gRodon2"
output <- if (length(args) >= 2) args[[2]] else "tests/fixtures/codon_reference.json"
source_files <- file.path(upstream, "R", c("utilities.R", "codon_pair_bias.R", "codon_usage_bias.R", "predict_growth.R"))
for (path in source_files) source(path)
load(file.path(upstream, "R", "sysdata.rda"))

set.seed(20261019)
codons <- sort(names(getGeneticCode("11")))
# Sample codons rather than organisms: these are mathematical test strings and
# may contain internal stop codons, intentionally exercising their treatment.
sequences <- vapply(seq_len(24), function(i) {
  probability <- if (i <= 8) seq_along(codons)^1.4 else rev(seq_along(codons))^0.7
  paste(sample(codons, 90 + 3 * i, replace = TRUE, prob = probability), collapse = "")
}, character(1))
he <- seq_along(sequences) <= 8
depth <- rep(c(0, 1, 2, 5, 13, 3), length.out = length(sequences))
metrics <- c("MILC", "ENCprime", "B", "SCUO", "MCB")

make_case <- function(name, seqs = sequences, high = he, options = list()) {
  input <- DNAStringSet(seqs)
  normal <- suppressWarnings(do.call(getCodonStatistics, c(list(genes = input, highly_expressed = high), options)))
  all <- suppressWarnings(do.call(getCodonStatistics, c(list(genes = input, highly_expressed = high, all_metrics = TRUE), options)))
  list(name = name, genes = unname(seqs), highly_expressed = unname(high),
       options = options, statistics = as.list(normal[1, ]), all_metrics = as.list(all[1, ]))
}

cases <- list()
for (code in c("1", "4", "11")) {
  cases[[length(cases) + 1]] <- make_case(paste0("genetic_code_", code), options = list(genetic_code = code))
  cases[[length(cases) + 1]] <- make_case(paste0("weighted_code_", code), options = list(genetic_code = code, depth_of_coverage = depth))
}
cases[[length(cases) + 1]] <- make_case("filtered_sequences", seqs = c(sequences, strrep("AAA", 80), strrep("AAN", 90), strrep("AAA", 90), "AAA"), high = c(he, FALSE, FALSE, TRUE, FALSE))
cases[[length(cases) + 1]] <- make_case("fragment_threshold", seqs = substring(sequences, 1, 150), options = list(fragments = TRUE))
cases[[length(cases) + 1]] <- make_case("trimmed_start", options = list(fragments = NA, trimlen = 150, trimside = "start"))
cases[[length(cases) + 1]] <- make_case("trimmed_with_short_retained", seqs = c(substring(sequences[1:4], 1, 90), sequences[5:24]), options = list(fragments = NA, trimlen = 150, trimside = "start"))
cases[[length(cases) + 1]] <- make_case("trimmed_short_excluded", seqs = c(substring(sequences[1:4], 1, 90), sequences[5:24]), options = list(fragments = NA, trimlen = 150, trimside = "start", exclude_short = TRUE))
cases[[length(cases) + 1]] <- make_case("weighted_trimmed", options = list(fragments = NA, trimlen = 240, trimside = "start", depth_of_coverage = depth))
restricted <- c(strrep("AAAAAG", 50), paste0(strrep("AAA", 80), strrep("AAG", 20)), paste0(strrep("AAA", 20), strrep("AAG", 80)))
cases[[length(cases) + 1]] <- make_case("restricted_amino_acids", seqs = restricted, high = c(TRUE, TRUE, FALSE))
cases[[length(cases) + 1]] <- make_case("balanced_all_codons", seqs = rep(paste(rep(codons, 10), collapse = ""), 3), high = c(TRUE, TRUE, FALSE))

# Low-level coRdon comparisons, including non-default translation conventions.
low_level <- list()
table <- codonTable(DNAStringSet(sequences[1:6]))
for (code in c("1", "4", "11")) {
  for (alternative_starts in c(TRUE, FALSE)) {
    expected <- lapply(metrics, function(method) {
      values <- do.call(get(method), list(cTobject = table, id_or_name2 = code, alt.init = alternative_starts))
      if (method == "SCUO") unname(values) else unname(values[, 1])
    })
    names(expected) <- metrics
    low_level[[length(low_level) + 1]] <- list(
      genetic_code = code, alternative_starts = alternative_starts,
      genes = unname(sequences[1:6]), expected = expected
    )
  }
}

# Store the actual pooled random backgrounds: Python's and R's RNG streams are
# intentionally different, but identical count inputs must yield identical metrics.
fixed_backgrounds <- lapply(seq_len(4), function(i) {
  gene <- sequences[[i]]
  bases <- strsplit(gene, "", fixed = TRUE)[[1]]
  shuffles <- replicate(100, paste(sample(bases, replace = FALSE), collapse = ""))
  tab <- codonTable(DNAStringSet(c(gene, shuffles)))
  expected <- lapply(metrics, function(method) {
    values <- do.call(get(method), list(cTobject = tab)) # coRdon default code 1
    if (method == "SCUO") unname(values[[1]]) else unname(values[1, 1])
  })
  names(expected) <- metrics
  list(gene = gene, reference_counts = unname(colSums(codonCounts(tab))), expected = expected)
})

pair_cases <- lapply(list(c("AAAAAG", "AAAAAG"), "TTGATGCTGGTG", sequences[1:3]), function(seqs) {
  list(genes = unname(seqs), genetic_code = "11", expected = getCPB(getPairCounts(DNAStringSet(seqs), genetic_code = "11")))
})

# Direct end-to-end predictGrowth calls test statistics, model dispatch and
# inverse transforms together. These inputs have ConsistencyHE < 0.6, so both
# metagenome_v2 paths stay deterministic without cross-language RNG matching.
stopifnot(getCodonStatistics(DNAStringSet(sequences), he)$ConsistencyHE < 0.6)
predictions <- list()
for (mode in c("full", "partial", "metagenome_v1", "metagenome_v2", "eukaryote")) {
  for (temperature in list("none", 27)) {
    options <- list(mode = mode, temperature = temperature, training_set = "madin")
    value <- suppressWarnings(do.call(predictGrowth, c(list(genes = DNAStringSet(sequences), highly_expressed = he), options)))
    predictions[[length(predictions) + 1]] <- list(name = paste(mode, temperature, sep = "_"), genes = unname(sequences), highly_expressed = unname(he), options = options, expected = value)
    if (mode %in% c("metagenome_v1", "metagenome_v2")) {
      options$depth_of_coverage <- depth
      value <- suppressWarnings(do.call(predictGrowth, c(list(genes = DNAStringSet(sequences), highly_expressed = he), options)))
      predictions[[length(predictions) + 1]] <- list(name = paste(mode, temperature, "weighted", sep = "_"), genes = unname(sequences), highly_expressed = unname(he), options = options, expected = value)
    }
  }
}

packages <- c("coRdon", "Biostrings", "matrixStats", "dplyr", "jsonlite")
versions <- lapply(packages, function(package) as.character(packageVersion(package)))
names(versions) <- packages
source_md5 <- as.list(unname(tools::md5sum(source_files)))
names(source_md5) <- basename(source_files)
fixture <- list(
  provenance = list(generator = "scripts/generate_codon_reference.R", R = R.version.string,
                    packages = versions, upstream_version = read.dcf(file.path(upstream, "DESCRIPTION"))[1, "Version"],
                    source_md5 = source_md5, model_data_md5 = unname(tools::md5sum(file.path(upstream, "R", "sysdata.rda"))), seed = 20261019L,
                    note = "Original gRodon R functions and installed coRdon executed on artificial sequences."),
  cases = cases, low_level = low_level, fixed_backgrounds = fixed_backgrounds,
  pair_cases = pair_cases, predictions = predictions,
  weighted_medians = list(
    list(values = c(0, 10), weights = c(1, 3), expected = weightedMedian(c(0, 10), c(1, 3))),
    list(values = c(10, 2, 5, 1, 3), weights = c(1, 4, 2, 0, 3), expected = weightedMedian(c(10, 2, 5, 1, 3), c(1, 4, 2, 0, 3)))
  )
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write_json(fixture, output, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null", null = "null")
cat("Wrote", length(cases), "aggregate cases,", length(low_level), "metric cases,", length(fixed_backgrounds), "fixed shuffle backgrounds and", length(predictions), "end-to-end predictions to", output, "\n")
