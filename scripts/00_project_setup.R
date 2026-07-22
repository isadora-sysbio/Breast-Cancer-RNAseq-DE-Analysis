# Project 2B: shared setup and safety helpers
# Run analysis scripts from anywhere inside the project directory.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    has_metadata <- file.exists(file.path(current, "metadata", "sample_sheet.csv"))
    has_results <- dir.exists(file.path(current, "results", "salmon"))

    if (has_metadata && has_results) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not find the Project 2B root. Expected metadata/sample_sheet.csv ",
        "and results/salmon/. Start R inside Breast-Cancer-RNAseq-DE-Analysis."
      )
    }
    current <- parent
  }
}

PROJECT_DIR <- find_project_root()
setwd(PROJECT_DIR)

ALLOW_OVERWRITE <- identical(
  toupper(Sys.getenv("PROJECT2B_ALLOW_OVERWRITE", unset = "NO")),
  "YES"
)

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Install them in the rnaseq environment before continuing."
    )
  }
}

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  invisible(path)
}

assert_dir <- function(path) {
  if (!dir.exists(path)) {
    stop("Required directory not found: ", path)
  }
  invisible(path)
}

prepare_output <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !ALLOW_OVERWRITE) {
    stop(
      "Output already exists and was protected from overwriting: ", path, "\n",
      "Run in a clean project copy, or only after verifying your backup set:\n",
      "export PROJECT2B_ALLOW_OVERWRITE=YES"
    )
  }
  invisible(path)
}

save_rds_safe <- function(object, path) {
  prepare_output(path)
  saveRDS(object, path)
}

write_csv_safe <- function(object, path, row.names = FALSE) {
  prepare_output(path)
  write.csv(object, path, row.names = row.names)
}

write_lines_safe <- function(text, path) {
  prepare_output(path)
  writeLines(text, path)
}

save_plot_safe <- function(plot, path, width = 10, height = 7, dpi = 300) {
  prepare_output(path)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

choose_column <- function(data, candidates, label, required = TRUE) {
  available <- names(data)
  exact <- candidates[candidates %in% available]

  if (length(exact) > 0) {
    return(exact[[1]])
  }

  lower_available <- tolower(available)
  lower_candidates <- tolower(candidates)
  matched_index <- match(lower_candidates, lower_available, nomatch = 0L)
  matched_index <- matched_index[matched_index > 0L]

  if (length(matched_index) > 0) {
    return(available[[matched_index[[1]]]])
  }

  if (required) {
    stop(
      "Could not identify the ", label, " column in metadata/sample_sheet.csv. ",
      "Expected one of: ", paste(candidates, collapse = ", "), ". ",
      "Available columns: ", paste(available, collapse = ", ")
    )
  }

  NULL
}

normalize_group <- function(x) {
  cleaned <- trimws(as.character(x))
  upper <- toupper(gsub("[- +]", "_", cleaned))

  result <- rep(NA_character_, length(cleaned))
  result[upper %in% c("ER_POS", "ER_POSITIVE", "ER+", "ERPOS")] <- "ER_pos"
  result[upper %in% c("TNBC", "TRIPLE_NEGATIVE", "TRIPLE_NEGATIVE_BREAST_CANCER")] <- "TNBC"

  result
}

load_sample_design <- function(path = "metadata/sample_sheet.csv") {
  assert_file(path)
  metadata <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  srr_column <- choose_column(
    metadata,
    c("SRR", "Run", "run", "run_accession", "SRA_Run", "accession"),
    "SRR/run accession"
  )
  group_column <- choose_column(
    metadata,
    c("Group", "group", "Condition", "condition", "Subtype", "subtype"),
    "biological group"
  )
  sample_column <- choose_column(
    metadata,
    c("Sample", "sample", "sample_id", "Sample_ID", "sample_name", "SampleName"),
    "sample ID",
    required = FALSE
  )

  srr <- trimws(as.character(metadata[[srr_column]]))
  sample_id <- if (is.null(sample_column)) {
    srr
  } else {
    trimws(as.character(metadata[[sample_column]]))
  }
  group <- normalize_group(metadata[[group_column]])

  if (anyNA(group)) {
    bad_values <- unique(metadata[[group_column]][is.na(group)])
    stop(
      "Unrecognized group value(s) in sample sheet: ",
      paste(bad_values, collapse = ", "),
      ". Expected ER_pos/ER+ or TNBC."
    )
  }
  if (any(sample_id == "") || any(srr == "")) {
    stop("Blank sample IDs or SRR accessions were found in the sample sheet.")
  }
  if (anyDuplicated(sample_id)) {
    stop("Sample IDs must be unique in metadata/sample_sheet.csv.")
  }
  if (anyDuplicated(srr)) {
    stop("SRR accessions must be unique in metadata/sample_sheet.csv.")
  }

  design <- data.frame(
    sample_id = sample_id,
    SRR = srr,
    Group = factor(group, levels = c("TNBC", "ER_pos")),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )

  if (nrow(design) != 20L) {
    warning("Expected 20 samples for Project 2B, but found ", nrow(design), ".")
  }
  if (!all(c("TNBC", "ER_pos") %in% unique(as.character(design$Group)))) {
    stop("Both TNBC and ER_pos samples are required for the planned comparison.")
  }

  design
}

for (directory in c(
  "counts", "docs", "figures", "logs", "results/deseq2",
  "results/pathway_inputs", "results/pathway_outputs",
  "results/reactome_outputs", "results/qc", "results/r_objects"
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

message("Project root: ", PROJECT_DIR)
message("Overwrite protection: ", if (ALLOW_OVERWRITE) "OFF" else "ON")
