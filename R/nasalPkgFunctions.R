########################################################
## Code: Marcelo Navarro Diaz
## Contact: marcelo.n.d@ciencias.unam.mx
########################################################
# ----- Functions for reading and handling data -----
#' Read metadata table from CSV
#'
#' Loads a metadata table from a CSV file, optionally sorting it by row names.
#' The first column of the CSV file is treated as row names (typically sample IDs).
#'
#' @param path Character. Path to the metadata CSV file.
#' @param sort_table Logical. If TRUE, metadata rows are sorted alphabetically
#'   by their row names (default = FALSE).
#'
#' @return A data.frame containing the metadata, with row names assigned from
#'   the first column of the CSV file.
#'
#' @details
#' The function expects the CSV file to have:
#' - A header row
#' - A first column containing unique sample identifiers
#'
#' @examples
#' \dontrun{
#' # Read metadata as-is
#' md <- read_metadata("metadata.csv")
#'
#' # Read metadata and sort by sample names
#' md_sorted <- read_metadata("metadata.csv", sort_table = TRUE)
#' }
#' @export
read_metadata <- function(path, sort_table = FALSE) {

  # Input checks
  if (!is.character(path) || length(path) != 1) {
    stop("`path` must be a single character string.")
  }
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  # Read CSV
  md <- read.csv(path, row.names = 1, stringsAsFactors = FALSE)

  # Optional: sort by row names
  if (isTRUE(sort_table)) {
    md <- md[order(row.names(md)), , drop = FALSE]
  }

  return(md)
}


#' Read feature table exported from FBMN/Hitchhiker's Guide workflow
#'
#' Imports a feature table generated following the workflow described in
#' Pakkir Shah et al. (Nature Protocols) for feature-based molecular networking.
#' The function reads the table, optionally sorts rows by feature names, removes
#' trailing `.mzML` from row names, and returns the **transposed** table so that
#' features become columns and samples become rows.
#'
#' @param path Character. Path to the exported feature table (CSV or delimited).
#' @param sort_by_names Logical. If TRUE, row names (usually sample names or
#'   feature IDs) are sorted alphabetically (default = FALSE).
#' @param p_sep Character. Field separator used in the file. Default is `","`.
#'
#' @return A transposed numeric data.frame where:
#'   - rows correspond to samples
#'   - columns correspond to features
#'   Feature IDs have `.mzML` removed.
#'
#' @details
#' The function expects an export format compatible with the Hitchhiker's Guide /
#' FBMN workflow (GNPS), typically including:
#' - a header row
#' - first column containing feature identifiers
#' - sample intensities across columns
#'
#' The returned object is transposed because many downstream analyses in this
#' project expect **samples as rows** and **features as columns**.
#'
#' @examples
#' \dontrun{
#' # Read feature table as exported by the FBMN workflow
#' ft <- read_ft("feature_table.csv")
#'
#' # Read and sort by row names
#' ft_sorted <- read_ft("feature_table.csv", sort_by_names = TRUE)
#' }
#' @export
read_ft <- function(path, sort_by_names = FALSE, p_sep = ",") {

  # ---- Input checks ----
  if (!is.character(path) || length(path) != 1) {
    stop("`path` must be a single character string.")
  }
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  if (!is.character(p_sep) || length(p_sep) != 1) {
    stop("`p_sep` must be a single character separator string.")
  }

  # Read feature table
  ft <- read.csv2(
    file   = path,
    header = TRUE,
    row.names = 1,
    sep    = p_sep,
    dec    = ".",
    stringsAsFactors = FALSE
  )

  # Optional: sort by row names
  if (isTRUE(sort_by_names)) {
    ft <- ft[order(row.names(ft)), , drop = FALSE]
  }

  # Clean row names (remove trailing .mzML)
  rownames(ft) <- gsub("\\.mzML$", "", rownames(ft))

  # Return transposed table
  # Samples become rows, features become columns.
  return(t(ft))
}

#' Filter features by minimum count across samples
#'
#' Filters rows (features) of a count table, keeping only those that have
#' at least `min_count` in at least `col_number` columns.
#'
#' @param feature_table A numeric matrix or data.frame with features in rows
#'   and samples in columns.
#' @param min_count Numeric scalar. Minimum count value for a column to be
#'   considered "present".
#' @param col_number Integer scalar. Minimum number of columns that must meet
#'   `min_count` for a feature to be kept.
#'
#' @return A subset of `feature_table` containing only rows that meet the
#'   filtering criteria. If no features pass the filter, an object with zero
#'   rows (but same columns) is returned.
#'
#' @examples
#' \dontrun{
#' # Keep features present (>= 10 counts) in at least 3 samples
#' filtered <- filter_features_by_col_counts(ft, min_count = 10, col_number = 3)
#' }
#' @export
filter_features_by_col_counts <- function(feature_table, min_count, col_number) {

  # Input checks
  if (!is.matrix(feature_table) && !is.data.frame(feature_table)) {
    stop("`feature_table` must be a matrix or data.frame.")
  }
  if (!is.numeric(min_count) || length(min_count) != 1) {
    stop("`min_count` must be a single numeric value.")
  }
  if (!is.numeric(col_number) || length(col_number) != 1) {
    stop("`col_number` must be a single numeric value.")
  }

  # Handle "no columns" edge case explicitly
  if (ncol(feature_table) == 0) {
    stop("`feature_table` has no columns.")
  }

  # Main filtering logic
  if (ncol(feature_table) > 1) {
    # Count how many columns per row are >= min_count, then filter
    keep <- rowSums(feature_table >= min_count) >= col_number
    return(feature_table[keep, , drop = FALSE])
  } else {
    # Single-column case: only keep rows where that column meets min_count
    ft <- feature_table[feature_table[, 1] >= min_count, , drop = FALSE]
    return(ft)
  }
}

#' Set selected species-sample combinations to zero
#'
#' Sets specified entries of a species-by-sample abundance matrix to zero.
#' This is useful for removing detected species from selected samples
#' (e.g., manual curation or experimental design adjustments).
#'
#' @param df A numeric matrix or data.frame where rows are species
#'   and columns are samples.
#' @param species_name Character. Single species name (must match a row name).
#' @param sample_names Character vector of sample names (must match column names).
#'
#' @return The same object as `df`, with the selected species in the selected
#'   samples set to zero.
#'
#' @details
#' The function performs two safety checks:
#'   - Species must exist in `rownames(df)`.
#'   - All provided sample names must exist in `colnames(df)`.
#'
#' If any name is missing, the function stops with an informative error.
#'
#' @examples
#' \dontrun{
#' # Remove species "Staphylococcus aureus" from samples SC7_1 and SC7_2
#' cleaned_df <- zero_out_species_in_samples(
#'   df = abundance_table,
#'   species_name = "Staphylococcus aureus",
#'   sample_names = c("SC7_1", "SC7_2")
#' )
#' }
#' @export
zero_out_species_in_samples <- function(df, species_name, sample_names) {

  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }
  if (!is.character(species_name) || length(species_name) != 1) {
    stop("`species_name` must be a single character string.")
  }
  if (!is.character(sample_names)) {
    stop("`sample_names` must be a character vector of column names.")
  }

  # Safety check: species must exist
  if (!(species_name %in% rownames(df))) {
    stop("Species '", species_name, "' not found in rownames(df).")
  }

  # Safety check: samples must exist
  missing_samples <- sample_names[!sample_names %in% colnames(df)]
  if (length(missing_samples) > 0) {
    stop("Samples not found in df: ",
         paste(missing_samples, collapse = ", "))
  }

  # Zero out selected species in selected samples
  df[species_name, sample_names] <- 0

  return(df)
}

#' Remove features whose names begin with specified prefixes
#'
#' Filters a feature table by removing any rows (features) whose row names
#' start with one or more given prefixes. Matching is performed using a
#' combined regular expression and is case-sensitive unless patterns
#' include case-insensitive constructs.
#'
#' @param df A matrix or data.frame with features in rows.
#' @param patterns Character vector of prefixes that should trigger removal.
#'   Each element is interpreted as a *literal prefix* used to identify
#'   unwanted features.
#'
#' @return The input `df` with all rows removed whose row names begin with
#'   any of the specified prefixes.
#'
#' @details
#' The function builds a single regex of the form:
#' \code{"^(prefix1|prefix2|...)"}
#' and removes all rows whose row names match this expression.
#'
#' If no rows match any prefix, the original table is returned unchanged.
#'
#' @examples
#' \dontrun{
#' # Remove all features beginning with "Blank" or "Control"
#' ft_clean <- remove_feature_by_prefix(
#'   df = feature_table,
#'   patterns = c("Blank", "Control")
#' )
#' }
#' @export
remove_feature_by_prefix <- function(df, patterns) {

  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("`patterns` must be a non-empty character vector.")
  }

  # Create a single regex pattern that matches any of the species names at the start
  combined_pattern <- paste0("^(", paste(patterns, collapse = "|"), ")")

  # Filter the dataframe: keep rows that do NOT match the pattern
  df_filtered <- df[!grepl(combined_pattern, rownames(df)), ]

  return(df_filtered)
}

#' Map species-level abundances to strain-level abundances
#'
#' Converts a species-level abundance table into a strain-level table using an
#' inoculation design matrix. For each sample, the abundance of a species in
#' `df1` is assigned to all strains of that species that were inoculated in the
#' corresponding SynCom (as specified in `df2`).
#'
#' @param df1 A numeric matrix or data.frame with **species-level abundances**.
#'   Rows are species, columns are samples (e.g., "SC7_1", "SC7_2", ...).
#' @param df2 A data.frame describing **strain inoculation**. The first column
#'   contains strain names (e.g., "Corynebacterium propinquum 16"), and the
#'   remaining columns are SynCom IDs (e.g., "SC7", "SC12", ...) with values
#'   0/1 indicating whether that strain was inoculated in that SynCom.
#'
#' @return A data.frame of strain-level abundances with:
#'   - rows = strains (from the first column of `df2`)
#'   - columns = samples (same as `df1` column names)
#'
#' @details
#' Species-strain mapping:
#' - Species names in `df1` are expected to match the **first two words** of the
#'   strain names in `df2` (e.g., `"Corynebacterium propinquum"` in `df1`
#'   matches `"Corynebacterium propinquum 16"` in `df2`).
#' - Species names are extracted from strain names using the pattern
#'   `^([A-Za-z]+ [A-Za-z]+).*`.
#'
#' Samples and SynComs:
#' - Sample names in `df1` (columns) are assumed to follow the pattern
#'   `"SCx_rep"` (e.g., `"SC7_1"`, `"SC7_2"`).
#' - The SynCom ID used to query `df2` is obtained by taking the part before
#'   the first underscore (e.g., `"SC7"` from `"SC7_1"`).
#'
#' Strains whose species are not present in `df1` remain at zero.
#'
#' @examples
#' \dontrun{
#' # df1: species-level OTU/abundance table
#' # df2: inoculation design (first column = strains, other columns = SynComs)
#' strain_abundance <- merge_abundance_by_strain(df1 = species_abundance,
#'                                               df2 = inoculation_design)
#' }
#' @export
merge_abundance_by_strain <- function(df1, df2) {

  # Input checks
  if (!is.matrix(df1) && !is.data.frame(df1)) {
    stop("`df1` must be a matrix or data.frame with species-level abundances.")
  }
  if (!is.matrix(df2) && !is.data.frame(df2)) {
    stop("`df2` must be a matrix or data.frame with strain design information.")
  }

  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)

  if (is.null(rownames(df1))) {
    stop("`df1` must have species names as rownames.")
  }
  if (ncol(df2) < 2) {
    stop("`df2` must have at least two columns: strains and SynCom indicators.")
  }

  # Helper: get inoculated strains for a given SynCom
  get_inoculated_strains <- function(df2, sample_name) {
    if (!sample_name %in% colnames(df2)) {
      stop("SynCom column '", sample_name, "' not found in df2.")
    }
    sample_column <- df2[[sample_name]]
    inoculated_indices <- which(sample_column == 1)
    df2[inoculated_indices, 1]  # first column = strain names
  }

  # Extract names
  species_names_df1 <- rownames(df1)
  strain_names_df2  <- df2[, 1]

  # Prepare empty strain-level abundance matrix
  new_abundance_matrix <- matrix(
    0,
    nrow = nrow(df2),
    ncol = ncol(df1)
  )
  rownames(new_abundance_matrix) <- strain_names_df2
  colnames(new_abundance_matrix) <- colnames(df1)

  samples <- colnames(new_abundance_matrix)

  # Optional consistency check: required SynCom IDs in df2
  syncom_ids_from_samples <- vapply(strsplit(samples, "_"), `[`, character(1), 1)
  missing_syncoms <- setdiff(unique(syncom_ids_from_samples), colnames(df2))
  if (length(missing_syncoms) > 0) {
    stop("SynCom IDs missing in df2: ",
         paste(missing_syncoms, collapse = ", "))
  }

  # Main loop over samples
  for (i in seq_along(samples)) {
    sample_name <- samples[i]

    # Extract SynCom ID from sample name (e.g., "SC7" from "SC7_1")
    current_sc <- strsplit(sample_name, "_")[[1]][1]

    # Get strains inoculated in this SynCom
    inoc_strains_per_sample <- get_inoculated_strains(df2 = df2,
                                                      sample_name = current_sc)

    # Loop over strains inoculated in this SynCom
    for (x in seq_along(inoc_strains_per_sample)) {
      strain_name <- inoc_strains_per_sample[x]

      # Row index in df2 / new_abundance_matrix
      index_strain_df2 <- which(strain_names_df2 == strain_name)

      # Derive species name from strain name (first two words)
      species_name <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", strain_name)

      if (species_name %in% species_names_df1) {
        index_species_df1 <- which(species_names_df1 == species_name)

        # Abundance of this species in this sample
        current_abundance <- df1[index_species_df1, sample_name]

        # Assign abundance to the corresponding strain row + sample column
        new_abundance_matrix[index_strain_df2, i] <- current_abundance
      }
    }
  }

  return(as.data.frame(new_abundance_matrix))
}

#' Merge non-target strains to species-level
#'
#' Takes a strain-level feature (OTU) table and collapses all **non-target**
#' strains to a single pseudo-strain per species, while keeping **target**
#' strains unchanged.
#'
#' @param df A numeric matrix or data.frame with strains in rows and samples
#'   in columns. Row names are expected to contain species and strain
#'   information, where the **first two words** encode the species name
#'   (e.g., "Corynebacterium propinquum 16").
#' @param target_species Character vector of species names (e.g.,
#'   `"Corynebacterium propinquum"`) for which **individual strains should be kept**.
#'
#' @return A data.frame where:
#'   - Rows corresponding to `target_species` remain at strain-level (unchanged).
#'   - All other strains are **aggregated by species**, and each species is
#'     represented by a single row whose name is `"<Genus> <species> 1"`.
#'
#' @details
#' Species names are derived from the row names by taking the first two words:
#' \code{paste(x[1:2], collapse = " ")} after splitting by spaces.
#'
#' For non-target species:
#'   - All strain rows for a given species are summed across rows
#'   - The resulting aggregated row is named `"<Genus> <species> 1"` to preserve
#'     compatibility with functions that expect a strain-like naming convention.
#'
#' If there are no non-target strains, only the target strain rows are returned.
#'
#' @examples
#' \dontrun{
#' # Keep all strains of C. propinquum separate, merge all others by species
#' merged_df <- merge_non_target_strains(
#'   df = strain_table,
#'   target_species = c("Corynebacterium propinquum")
#' )
#' }
#' @export
merge_non_target_strains <- function(df, target_species) {

  # Input checks
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame with strains in rows.")
  }
  if (is.null(rownames(df))) {
    stop("`df` must have row names containing species and strain information.")
  }
  if (!is.character(target_species)) {
    stop("`target_species` must be a character vector of species names.")
  }

  df <- as.data.frame(df)

  # Extract species names from rownames (first two words)
  species_names <- sapply(
    strsplit(rownames(df), " "),
    function(x) paste(x[1:2], collapse = " ")
  )

  # Identify target vs non-target rows
  is_target <- species_names %in% target_species

  target_df     <- df[is_target, , drop = FALSE]
  non_target_df <- df[!is_target, , drop = FALSE]
  non_target_species <- species_names[!is_target]

  # Aggregate non-target strains by species
  if (nrow(non_target_df) > 0) {
    aggregated <- aggregate(
      non_target_df,
      by = list(Species = non_target_species),
      FUN = sum
    )

    # Row names become "<Genus> <species> 1"
    rownames(aggregated) <- paste(aggregated$Species, "1")
    aggregated$Species <- NULL
  } else {
    aggregated <- NULL
  }

  # Combine target strains with aggregated non-target species
  result <- rbind(target_df, aggregated)

  return(result)
}


#' Cluster samples based on relative abundance profiles
#'
#' Performs sample clustering on a feature (OTU) table using relative abundance,
#' hierarchical clustering (Ward.D2), and PAM clustering with silhouette-based
#' estimation of the optimal number of clusters (K) when not provided.
#'
#' @param abundance_df A numeric matrix or data.frame with features (e.g. OTUs,
#'   species) in rows and samples in columns.
#' @param k Optional integer. Number of clusters to use for PAM. If \code{NULL}
#'   (default), K is estimated using the average silhouette width over
#'   K = 2,...,K_max, where K_max = \code{min(10, n_samples - 1)}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{clusters}{A data.frame with columns \code{Sample} and \code{Cluster}
#'     containing the PAM cluster assignment for each sample.}
#'   \item{best_k}{The number of clusters used (either provided via \code{k} or
#'     estimated using the silhouette method).}
#'   \item{rel_abundance_ordered}{A data.frame of relative abundances with
#'     samples ordered according to the hierarchical clustering dendrogram.}
#'   \item{sample_order}{Character vector of sample names in dendrogram order.}
#' }
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Convert the input to a matrix.
#'   \item Convert counts to relative abundance per sample
#'         (each column sums to 1).
#'   \item Compute a distance matrix between samples using Euclidean distance
#'         on transposed relative abundance (\code{dist(t(mat_rel))}).
#'   \item Perform hierarchical clustering with method \code{"ward.D2"}.
#'   \item If \code{k} is \code{NULL}, estimate the best K by fitting PAM
#'         clustering models for K = 2,...,K_max and choosing the K that
#'         maximizes the average silhouette width.
#'   \item Fit the final PAM model with \code{best_k} clusters and extract the
#'         cluster assignments.
#'   \item Order samples according to the hierarchical clustering dendrogram and
#'         reorder the relative abundance matrix accordingly.
#' }
#'
#' Requires the \pkg{cluster} package for \code{cluster::pam()}.
#'
#' @examples
#' \dontrun{
#' res <- cluster_samples(abundance_df = otu_table)
#' }
#' @export
cluster_samples <- function(abundance_df, k = NULL) {

  # Input checks
  if (!is.matrix(abundance_df) && !is.data.frame(abundance_df)) {
    stop("`abundance_df` must be a numeric matrix or data.frame.")
  }

  mat <- as.matrix(abundance_df)

  if (!is.numeric(mat)) {
    stop("`abundance_df` must be numeric (counts or abundances).")
  }

  n_samples <- ncol(mat)
  if (n_samples < 2) {
    stop("`abundance_df` must contain at least 2 samples (columns) for clustering.")
  }

  # Relative abundance per sample
  col_sums <- colSums(mat, na.rm = TRUE)
  if (any(col_sums == 0)) {
    stop("One or more samples have total abundance of zero; cannot compute relative abundances.")
  }

  mat_rel <- sweep(mat, 2, col_sums, "/")

  mat_input <- mat_rel

  # Distance matrix and hierarchical clustering
  dist_mat <- stats::dist(t(mat_input), method = "euclidean")
  hc <- stats::hclust(dist_mat, method = "ward.D2")

  # Determine number of clusters (best_k)
  if (is.null(k)) {
    # Upper bound for K: at most 10 or n_samples - 1
    k_max <- min(10, n_samples - 1)

    if (k_max < 2) {
      stop("Not enough samples to estimate clusters (need at least 3 samples).")
    }

    ks_to_try <- 2:k_max

    sil_widths <- sapply(ks_to_try, function(k_try) {
      pam_fit <- cluster::pam(dist_mat, k_try)
      pam_fit$silinfo$avg.width
    })

    best_k <- ks_to_try[which.max(sil_widths)]
  } else {
    if (!is.numeric(k) || length(k) != 1 || k <= 1) {
      stop("`k` must be a single integer greater than 1.")
    }
    if (k >= n_samples) {
      stop("`k` must be less than the number of samples.")
    }
    best_k <- as.integer(k)
  }

  # Final PAM clustering with best_k
  pam_fit <- cluster::pam(dist_mat, best_k)
  clusters <- pam_fit$clustering

  cluster_df <- data.frame(
    Sample  = names(clusters),
    Cluster = as.factor(clusters),
    row.names = NULL
  )

  # Reorder samples by dendrogram order
  sample_order <- hc$labels[hc$order]
  mat_rel_ordered <- mat_rel[, sample_order, drop = FALSE]

  # Return results
  return(list(
    clusters = cluster_df,
    best_k = best_k,
    rel_abundance_ordered = as.data.frame(mat_rel_ordered),
    sample_order = sample_order
  ))
}

#' Transform a feature table using scaling or normalization methods
#'
#' Applies one of three common transformations to a feature table:
#' \describe{
#'   \item{\code{"zscale"}}{Z-score standardization of each column using \code{scale()}.}
#'   \item{\code{"min_max"}}{Min-max normalization of each numeric column to the range between 0 and 1.}
#'   \item{\code{"rel_abundance"}}{Conversion to relative abundance where each column sums to 1.}
#' }
#'
#' @param feature_table A numeric matrix or data.frame with features in rows
#'   and samples in columns.
#' @param transform_method Character string. One of:
#'   \code{"zscale"}, \code{"min_max"}, \code{"rel_abundance"}.
#'
#' @return A transformed data.frame with the same dimensions as \code{feature_table}.
#'
#' @details
#' \itemize{
#'   \item For \code{"zscale"}, columns are centered and scaled to unit variance.
#'   \item For \code{"min_max"}, non-numeric columns are left unchanged.
#'   \item For \code{"rel_abundance"}, columns with sum zero will trigger an error
#'         to avoid division by zero.
#' }
#'
#' @examples
#' \dontrun{
#' scaled  <- transform_feature_table(otu_table, "zscale")
#' normed  <- transform_feature_table(otu_table, "min_max")
#' relab   <- transform_feature_table(otu_table, "rel_abundance")
#' }
#' @export
transform_feature_table <- function(feature_table, transform_method) {

  # ---- Input checks ----
  if (!is.matrix(feature_table) && !is.data.frame(feature_table)) {
    stop("`feature_table` must be a matrix or data.frame.")
  }

  if (!transform_method %in% c("zscale", "min_max", "rel_abundance")) {
    stop("Invalid `transform_method`. Must be 'zscale', 'min_max', or 'rel_abundance'.")
  }

  feature_table <- as.data.frame(feature_table)

  # ---- Transformation logic ----
  if (transform_method == "zscale") {

    df_transformed <- as.data.frame(scale(feature_table))

  } else if (transform_method == "min_max") {

    df_transformed <- feature_table
    numeric_cols <- sapply(df_transformed, is.numeric)

    normalize <- function(x) {
      rng <- max(x) - min(x)
      if (rng == 0) {
        warning("A column has zero variance; min-max scaling will return zeros.")
        return(rep(0, length(x)))
      }
      (x - min(x)) / rng
    }

    df_transformed[numeric_cols] <- lapply(df_transformed[numeric_cols], normalize)

  } else if (transform_method == "rel_abundance") {

    col_sums <- colSums(feature_table, na.rm = TRUE)
    if (any(col_sums == 0)) {
      stop("One or more columns sum to 0; cannot compute relative abundance normalization.")
    }

    df_transformed <- sweep(feature_table, 2, col_sums, FUN = "/")
    df_transformed <- as.data.frame(df_transformed)
  }

  return(df_transformed)
}

#' Sort Nanopore feature table columns by barcode order
#'
#' Sorts the columns of a Nanopore feature/OTU table by barcode-like column
#' names, using a mixed strategy that orders first by the character length of
#' the column name and then lexicographically. This is useful when barcodes are
#' named like "BC1", "BC2", ..., "BC10", so that "BC2" comes before "BC10".
#'
#' Optionally, new column names can be assigned after sorting.
#'
#' @param df A matrix or data.frame with features in rows and Nanopore samples
#'   (barcodes) in columns.
#' @param new_names Optional character vector of new column names. If provided,
#'   its length must match the number of columns in \code{df}.
#'
#' @return A data.frame with columns sorted by barcode order. If
#'   \code{new_names} is provided, columns are renamed accordingly.
#'
#' @details
#' Column names are sorted using:
#' \itemize{
#'   \item primary key: \code{nchar(colname)}
#'   \item secondary key: lexicographic order of \code{colname}
#' }
#'
#' This avoids the usual problem where purely lexicographic sorting would place
#' "BC10" before "BC2".
#'
#' @examples
#' \dontrun{
#' # Sort a Nanopore OTU table by barcode-like column names
#' df_sorted <- sort_nanopore_table_by_barcodes(nanopore_otu)
#'
#' # Sort and assign human-readable sample names
#' df_sorted_named <- sort_nanopore_table_by_barcodes(
#'   df = nanopore_otu,
#'   new_names = paste0("Sample_", seq_len(ncol(nanopore_otu)))
#' )
#' }
#' @export
sort_nanopore_table_by_barcodes <- function(df, new_names = NULL) {

  # ---- Input checks ----
  if (!is.matrix(df) && !is.data.frame(df)) {
    stop("`df` must be a matrix or data.frame.")
  }

  df <- as.data.frame(df)

  if (is.null(colnames(df))) {
    stop("`df` must have column names (barcodes).")
  }

  if (!is.null(new_names)) {
    if (!is.character(new_names)) {
      stop("`new_names` must be a character vector.")
    }
    if (length(new_names) != ncol(df)) {
      stop("Length of `new_names` (", length(new_names),
           ") must match the number of columns in `df` (", ncol(df), ").")
    }
  }

  # ---- Sort column names by length, then lexicographically ----
  cn <- colnames(df)
  sorted_names <- cn[order(nchar(cn), cn)]

  df_sorted <- df[, sorted_names, drop = FALSE]

  # ---- Optionally rename columns ----
  if (!is.null(new_names)) {
    colnames(df_sorted) <- new_names
  }

  return(df_sorted)
}

# ----- Cluster Barplots -----
#' Generate a Palette of Distinct Colors
#'
#' Returns a vector of distinct colors for plotting. By default, the function
#' samples \code{nColors} colors from a predefined palette of visually distinct
#' color values. Colors may be sampled with or without replacement.
#'
#' @param nColors Integer. Number of colors to return. Default is \code{60}.
#' @param replace_cols Logical. Whether to sample colors with replacement.
#'   Default is \code{FALSE}.
#'
#' @return A character vector of color hex codes or named R colors of length
#'   \code{nColors}.
#'
#' @examples
#' # Generate 10 distinct colors
#' get_palette(10)
#'
#' # Sample colors with replacement
#' get_palette(10, replace_cols = TRUE)
#'
#' @export
get_palette <- function(nColors = 60, replace_cols = FALSE){
  colors_vec <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3",
    "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
    "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
    "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
    "lavenderblush4", "lawngreen", "indianred1", "lightblue1", "honeydew4",
    "hotpink", "#e3ae78", "#a23f3f", "#290f76", "#ce7e00",
    "#386857", "#738564", "#e89d56", "#cd541d", "#1a3a46",
    "#9C4A1A", "#ffe599", "#583E26", "#A78B71", "#F7C815",
    "#EC9704", "#4B1E19", "firebrick2", "#C8D2D1", "#14471E",
    "#6279B8", "#DA6A00", "#C0587E", "#FC8B5E", "#FEF4C0",
    "#EA592A", "khaki3", "lavenderblush3", "indianred4", "lightblue",
    "honeydew1", "hotpink4", "ivory3", "#49516F", "#502F4C",
    "#A8C686", "#669BBC", "#29335C", "#E4572E", "#F3A712",
    "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD"
  )

  sample(colors_vec, nColors, replace = replace_cols)
}


#' Create Relative Abundance Barplots by Cluster
#'
#' Creates a stacked relative abundance barplot with one panel per cluster.
#' Samples are grouped by their assigned cluster and displayed as stacked bars
#' of bacterial relative abundance. Optionally, strain-level information can be
#' encoded using fill colors for species and patterns for strains.
#'
#' @param abundance_df A matrix or data frame of relative abundances with
#'   taxa (e.g., species or strain labels) in rows and samples in columns.
#'   Row names are used as the \code{Bacteria} identifiers in the plot.
#' @param cluster_df A data frame containing at least two columns:
#'   \code{Sample} (sample IDs matching the column names of \code{abundance_df})
#'   and \code{Cluster} (cluster assignment for each sample).
#' @param sample_order Optional character vector specifying the order of
#'   samples (column names of \code{abundance_df}) to be used on the x-axis.
#'   If \code{NULL}, the original column order of \code{abundance_df} is used.
#' @param colour_palette Optional named character vector of colors used to
#'   fill taxa. Names must match the \code{Bacteria} values (i.e., row names of
#'   \code{abundance_df}). If \code{NULL}, ggplot2's default color palette is
#'   used.
#' @param strains Logical. If \code{TRUE}, the function assumes that
#'   \code{Bacteria} names encode strain information in the last numeric token
#'   (e.g., \code{"Corynebacterium propinquum 1"}) and will represent species
#'   by fill color and strains by bar patterns. Default is \code{FALSE}.
#' @param best_k Integer or character. The number of clusters (k) to display in
#'   the plot title. This argument must be provided; otherwise, the function
#'   will stop with an error.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Orders the abundance matrix according to \code{sample_order}
#'     if provided.
#'   \item Optionally converts strain names to numeric strain IDs using
#'     \code{strain_name2strain_number()} when \code{strains = TRUE}.
#'   \item Reshapes the abundance matrix to long format and merges it with
#'     \code{cluster_df} via the \code{Sample} column.
#'   \item When \code{strains = FALSE}, creates a standard stacked barplot of
#'     taxa relative abundances per sample, faceted by cluster.
#'   \item When \code{strains = TRUE}, uses \pkg{ggpattern} to encode species
#'     as fill colors and strains as patterns within each stacked bar.
#' }
#'
#' The function relies on \pkg{ggplot2}, \pkg{reshape2}, \pkg{dplyr},
#' \pkg{ggpattern}, and a user-defined helper
#' \code{strain_name2strain_number()} that converts strain names to numeric
#' strain labels.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} object representing the relative abundance
#'     barplot faceted by cluster.}
#'   \item{df_long}{A data frame in long format containing the abundance data,
#'     cluster annotations, and (optionally) species/strain columns used for
#'     plotting.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage (no strain patterns, default colors)
#' res <- cluster_barplot_panels(
#'   abundance_df = abundance_mat,
#'   cluster_df   = sample_clusters,
#'   best_k       = 4
#' )
#' res$plot
#'
#' # With custom sample order and custom colors
#' res <- cluster_barplot_panels(
#'   abundance_df   = abundance_mat,
#'   cluster_df     = sample_clusters,
#'   sample_order   = c("S1", "S2", "S3"),
#'   colour_palette = my_colors,
#'   best_k         = 4
#' )
#'
#' # With strain-level encoding (requires strain_name2strain_number())
#' res <- cluster_barplot_panels(
#'   abundance_df = abundance_mat,
#'   cluster_df   = sample_clusters,
#'   strains      = TRUE,
#'   best_k       = 4
#' )
#' }
#'
#' @export
cluster_barplot_panels <- function(
    abundance_df,
    cluster_df,
    sample_order   = NULL,
    colour_palette = NULL,
    strains        = FALSE,
    best_k         = NULL
) {
  #require(cluster)
  #require(ggplot2)
  #require(reshape2)
  #require(dplyr)
  #require(ggpattern)

  if (is.null(best_k)) {
    stop("Argument 'best_k' must be provided.")
  }

  # Base matrix (all samples), optionally reordered
  mat_rel_ordered <- as.matrix(abundance_df)
  if (!is.null(sample_order)) {
    mat_rel_ordered <- mat_rel_ordered[, sample_order, drop = FALSE]
  }

  if (isTRUE(strains)) {
    message("Using strain data")
    # Convert table with strain names to a strain-number table
    mat_rel_ordered <- strain_name2strain_number(mat_rel_ordered)
  }

  # Melt for ggplot
  df_long <- reshape2::melt(mat_rel_ordered)
  colnames(df_long) <- c("Bacteria", "Sample", "Abundance")
  df_long <- merge(df_long, cluster_df, by = "Sample")

  # Add strain data columns if needed
  if (isTRUE(strains)) {
    df_long <- df_long %>%
      dplyr::mutate(
        strain   = paste0("Strain ", sub(".* ", "", Bacteria)),  # last token as strain
        species2 = sub(" \\d+$", "", Bacteria)                   # drop trailing number
      )
  }

  if (isFALSE(strains)) {
    p1 <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Bacteria)) +
      geom_bar(stat = "identity") +
      facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ylab("Relative Abundance") +
      ggtitle(paste("Stacked Barplot with Clusters (k =", best_k, ")"))
    message("Created plot without strain data")
  } else if (isTRUE(strains)) {
    # Clean the long-format table
    df_long <- df_long %>%
      dplyr::filter(!is.na(Abundance) & Abundance != 0) %>%
      dplyr::filter(!is.na(strain) & strain != 0)

    p1 <- ggplot(data = df_long, aes(x = Sample, y = Abundance)) +
      ggpattern::geom_bar_pattern(
        aes(fill = species2, pattern = strain),
        position        = "fill",
        stat            = "identity",
        show.legend     = TRUE,
        pattern_spacing = unit(2.5, "mm"),
        pattern_density = 0.0050,
        pattern_color   = "white",
        pattern_fill    = "white",
        pattern_angle   = 45
      ) +
      facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
      ggpattern::scale_pattern_manual(
        values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")
      ) +
      ggpattern::scale_pattern_spacing_manual(
        values = c(0, unit(0.025, "mm"), unit(0.025, "mm"))
      ) +
      guides(
        pattern = guide_legend(override.aes = list(fill = "grey")),
        fill    = guide_legend(override.aes = list(pattern = "none"))
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    message("Created plot with strain data")
  }

  if (!is.null(colour_palette)) {
    # Expecting a named vector: names must match Bacteria (rownames(abundance_df))
    p1 <- p1 + scale_fill_manual(values = colour_palette, drop = FALSE)
    message("Added custom color scale")
  }

  return(list(
    plot   = p1,
    df_long = df_long
  ))
}


#' Convert Strain-Level OTU/ASV Names to Species-Level Numbered Labels
#'
#' Converts row names of a data frame from full strain-level identifiers
#' (e.g., \code{"Genus species strainX"}) to a standardized species + strain ID
#' format (e.g., \code{"Genus species 1"}, \code{"Genus species 2"}).
#'
#' The function assumes that row names contain at least three whitespace-
#' separated components, where the first two correspond to the genus and species,
#' and the final component encodes the strain identity. Rows belonging to the
#' same species receive sequential numeric identifiers in the order they appear.
#'
#' @param df A data frame or matrix whose row names contain strain-level OTU/ASV
#'   identifiers. The row names must follow the structure:
#'   \code{"Genus species strainInfo"}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts the first two words (\code{"Genus species"}) from each row name.
#'   \item Assigns sequential numeric strain IDs within each species.
#'   \item Replaces the original row names with the standardized format
#'     \code{"Genus species <ID>"}.
#' }
#'
#' @return The input data frame with updated row names representing species-level
#'   groups with numeric strain identifiers.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(a = 1:3, b = 4:6)
#' rownames(df) <- c("Corynebacterium propinquum A1",
#'                   "Corynebacterium propinquum B7",
#'                   "Staphylococcus aureus T3")
#'
#' strain_name2strain_number(df)
#' }
#' @export
strain_name2strain_number <- function(df){
  # Extract only the "Genus species" part
  species_names <- sub(" \\S+$", "", rownames(df))

  # Create a numeric ID for each strain within the same species
  species_ids <- ave(species_names, species_names, FUN = function(x) seq_along(x))

  # Create new row names with species + strain ID
  new_rownames <- paste(species_names, species_ids)

  # Assign the new rownames to the dataframe
  rownames(df) <- new_rownames

  return(df)
}

#' Cluster Samples and Compute Mean Abundance of a Given Species
#'
#' Performs hierarchical clustering of samples based on their abundance profiles
#' and computes the mean relative abundance of a specified species within each
#' cluster. Optionally prints the sample names assigned to each cluster.
#'
#' @param df A numeric data frame or matrix where rows represent taxa (e.g.,
#'   species or strains) and columns represent samples.
#' @param species_name Character string giving the row name of the species whose
#'   mean abundance should be computed within each cluster.
#' @param k Integer. Number of clusters to cut the hierarchical clustering tree
#'   into. Default is \code{2}.
#' @param method Character string specifying the distance metric passed to
#'   \code{\link{dist}}. Default is \code{"euclidean"}.
#' @param show_samples Logical. If \code{TRUE}, prints the sample names belonging
#'   to each cluster. Default is \code{FALSE}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Computes a distance matrix between samples (\code{dist(t(df))}).
#'   \item Performs hierarchical clustering with \code{\link{hclust}}.
#'   \item Cuts the tree into \code{k} clusters using \code{\link{cutree}}.
#'   \item Computes the mean abundance of the target species within each cluster.
#' }
#'
#' Results (cluster means and optionally sample lists) are printed to the console.
#'
#' @return This function is called for its side effects (printed output) and does
#'   not return a value.
#'
#' @examples
#' \dontrun{
#' cluster_mean_abundance(df = abundance_matrix,
#'                        species_name = "Corynebacterium propinquum",
#'                        k = 3,
#'                        method = "euclidean",
#'                        show_samples = TRUE)
#' }
#'
#' @export
cluster_mean_abundance <- function(df, species_name, k = 2, method = "euclidean", show_samples = FALSE) {
  # Check species
  if (!(species_name %in% rownames(df))) {
    stop("Species not found in the dataframe.")
  }

  # Transpose for clustering samples
  dist_matrix <- dist(t(df), method = method)
  hc <- hclust(dist_matrix)

  # Cut tree into k groups
  groups <- cutree(hc, k = k)

  # Group sample names by cluster
  cluster_samples <- split(names(groups), groups)

  # Print results
  cat("Number of clusters:", length(cluster_samples), "\n\n")

  for (i in seq_along(cluster_samples)) {
    samples <- cluster_samples[[i]]
    mean_abund <- mean(as.numeric(df[species_name, samples]))

    cat("Cluster", i, "- Mean relative abundance of",
        species_name, ":", round(mean_abund, 5), "\n")

    if (show_samples) {
      cat("  Samples in cluster", i, ":\n")
      cat("   ", paste(samples, collapse = ", "), "\n\n")
    }
  }
}


#' Add a Cluster Column to a Metadata Data Frame
#'
#' Maps cluster assignments from a lookup table onto a metadata data frame by
#' matching key columns. The function performs sanity checks, warns about
#' duplicated keys in the cluster table, optionally warns about metadata rows
#' with missing matches, and appends a new column containing the mapped cluster
#' values.
#'
#' @param meta_df A data frame containing metadata. A key column within this
#'   data frame will be used to map cluster values onto it.
#' @param clusters_df A data frame containing cluster assignments. Must include
#'   a key column and a value column.
#' @param meta_key_col Character string. Name of the column in \code{meta_df}
#'   used to match entries against \code{clusters_df}.
#' @param cluster_key_col Character string. Name of the key column in
#'   \code{clusters_df} that corresponds to \code{meta_key_col}.
#' @param cluster_value_col Character string. Name of the column in
#'   \code{clusters_df} containing the cluster labels or values to be mapped.
#' @param new_col_name Character string. Name of the new column to be added to
#'   \code{meta_df} containing the mapped cluster entries.
#' @param warn_missing Logical. If \code{TRUE} (default), warns when some values
#'   in \code{meta_df[[meta_key_col]]} have no corresponding match in
#'   \code{clusters_df}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Checks that all required columns exist in both data frames.
#'   \item Warns if duplicate keys appear in the cluster table; in such cases,
#'         the last occurrence is used for mapping.
#'   \item Constructs a named lookup vector mapping keys to cluster values.
#'   \item Uses this lookup to assign cluster values to the metadata table.
#'   \item Optionally warns about unmatched keys in the metadata.
#' }
#'
#' This is a general-purpose utility for augmenting metadata with additional
#' annotations derived from separate data tables.
#'
#' @return A modified version of \code{meta_df} containing a new column
#'   \code{new_col_name} with mapped cluster values.
#'
#' @examples
#' \dontrun{
#' meta <- data.frame(Sample = c("S1", "S2", "S3"))
#' clusters <- data.frame(SampleID = c("S1", "S2"), Cluster = c(1, 2))
#'
#' add_cluster_column(
#'   meta_df          = meta,
#'   clusters_df      = clusters,
#'   meta_key_col     = "Sample",
#'   cluster_key_col  = "SampleID",
#'   cluster_value_col = "Cluster",
#'   new_col_name     = "Cluster_Assignment"
#' )
#' }
#'
#' @export
add_cluster_column <- function(meta_df,
                               clusters_df,
                               meta_key_col,
                               cluster_key_col,
                               cluster_value_col,
                               new_col_name,
                               warn_missing = TRUE) {
  # Basic checks
  for (nm in c(meta_key_col)) {
    if (!nm %in% names(meta_df)) stop(sprintf("Column '%s' not found in meta_df.", nm))
  }
  for (nm in c(cluster_key_col, cluster_value_col)) {
    if (!nm %in% names(clusters_df)) stop(sprintf("Column '%s' not found in clusters_df.", nm))
  }

  # Optional: detect duplicate keys in the cluster table
  dups <- duplicated(clusters_df[[cluster_key_col]])
  if (any(dups)) {
    dup_vals <- unique(clusters_df[[cluster_key_col]][dups])
    warning(sprintf(
      "Duplicate keys in clusters_df$%s for: %s. Using the last occurrence.",
      cluster_key_col, paste(head(dup_vals, 10), collapse = ", ")
    ))
  }

  # Build a named lookup vector: names = keys, values = clusters
  lookup <- setNames(clusters_df[[cluster_value_col]],
                     as.character(clusters_df[[cluster_key_col]]))

  # Map onto meta_df using the key column
  keys <- as.character(meta_df[[meta_key_col]])
  mapped <- unname(lookup[keys])

  # Warn if some keys were not found
  if (warn_missing) {
    missing_idx <- which(is.na(mapped) & !is.na(keys))
    if (length(missing_idx) > 0) {
      missing_vals <- unique(keys[missing_idx])
      warning(sprintf(
        "No cluster found for %d key(s) in meta_df$%s. Examples: %s",
        length(missing_vals), meta_key_col,
        paste(head(missing_vals, 10), collapse = ", ")
      ))
    }
  }

  # Add the new column
  meta_df[[new_col_name]] <- mapped
  meta_df
}

#' Order Samples According to Hierarchical Clustering
#'
#' Computes a hierarchical clustering of samples based on their feature
#' profiles and returns the sample names in the order determined by the
#' clustering dendrogram. This is useful for ordering heatmaps, barplots,
#' or other visualizations consistently with sample similarity.
#'
#' @param feature_table A numeric data frame or matrix where rows represent
#'   features (e.g., taxa, genes, metabolites) and columns represent samples.
#'   Row names are treated as feature identifiers.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Converts row names into a column and transposes the data so that
#'         samples become rows for clustering.
#'   \item Computes Euclidean distances between samples.
#'   \item Performs hierarchical clustering using Ward's \code{D2} method.
#'   \item Extracts sample names in the order they appear in the dendrogram.
#' }
#'
#' The resulting sample order is often used for heatmap column ordering or to
#' ensure consistent clustering-based visualization across analyses.
#'
#' @return A character vector of sample names ordered according to the
#'   hierarchical clustering.
#'
#' @examples
#' \dontrun{
#' sample_order <- order_samples_by_clustering(feature_table)
#' heatmap(feature_table[, sample_order])
#' }
#'
#' @export
order_samples_by_clustering <- function(feature_table){
  # Add feature names as a column (Species)
  df_otu <- feature_table %>% rownames_to_column(var = "Species")

  # Transpose table so samples become rows (remove Species column first)
  df_t <- as.matrix(t(df_otu[, -1]))

  # Hierarchical clustering of samples
  d <- dist(df_t, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")

  # Extract ordered sample names (skip the Species column name)
  ordered_samples_cluster <- colnames(df_otu)[-1][hc$order]

  return(ordered_samples_cluster)
}

#' Create a Grid of Relative Abundance Barplots Across Experiments
#'
#' Builds a faceted grid of stacked relative abundance barplots from multiple
#' feature tables (e.g., OTU/ASV count tables or relative abundance tables),
#' typically corresponding to different experiments or treatments.
#'
#' Each feature table is reshaped to long format, combined into a single data
#' frame, and plotted as stacked bars either:
#' \itemize{
#'   \item per sample, faceted by experiment (default), or
#'   \item per experiment, faceted by sample (when \code{shared_samples = TRUE}).
#' }
#'
#' Optionally, strain-level information can be encoded using patterns
#' (\pkg{ggpattern}) while species are encoded by fill color.
#'
#' @param feature_tables A list of numeric matrices or data frames, where rows
#'   represent features (e.g., taxa/strains) and columns represent samples. Each
#'   element of the list corresponds to one experiment or condition.
#' @param experiments_names Character vector of the same length as
#'   \code{feature_tables}, giving the experiment/treatment name for each
#'   feature table. These names are used to annotate and facet the plots.
#' @param shared_samples Logical. If \code{TRUE}, assumes that the same samples
#'   are shared across experiments and plots \code{experiment} on the x-axis
#'   with \code{sample} as facets. If \code{FALSE} (default), plots
#'   \code{sample} on the x-axis and facets by \code{experiment}.
#' @param strains Logical. If \code{TRUE}, treats row names of each feature
#'   table as strain-level identifiers of the form
#'   \code{"Genus species strainID"} and uses \code{\link{strain_name2strain_number}}
#'   and pattern aesthetics (\pkg{ggpattern}) to represent strains. Default
#'   is \code{FALSE}.
#' @param plot_title Character string. Overall title for the plot. Default is
#'   an empty string.
#' @param plot_title_size Numeric. Text size for the plot title. Default is
#'   \code{14}.
#' @param x_axis_text_size Numeric. Text size for x-axis tick labels. Default
#'   is \code{12}.
#' @param x_axis_title_size Numeric. Text size for the x-axis title. Default
#'   is \code{12}.
#' @param x_axis_text_angle Numeric. Rotation angle (in degrees) for x-axis
#'   tick labels. Default is \code{0}.
#' @param y_axis_title_size Numeric. Text size for the y-axis title. Default
#'   is \code{12}.
#' @param y_axis_text_size Numeric. Text size for y-axis tick labels. Default
#'   is \code{12}.
#' @param y_axis_text_angle Numeric. Rotation angle (in degrees) for y-axis
#'   tick labels. Default is \code{0}.
#' @param legend_pos Character string specifying the legend position (e.g.,
#'   \code{"right"}, \code{"bottom"}, \code{"none"}). Passed to
#'   \code{theme(legend.position = ...)}. Default is \code{"right"}.
#' @param legend_title_size Numeric. Text size for legend title. Default is
#'   \code{12}.
#' @param legend_text_size Numeric. Text size for legend text. Default is
#'   \code{12}.
#' @param legend_cols Integer. Number of columns in the legend for the
#'   \code{fill} aesthetic. Default is \code{3}.
#' @param legend_key_size Numeric. Size (in cm) of the legend key boxes.
#'   Default is \code{1}.
#' @param colour_palette Optional named character vector of colors used for
#'   species (or species groups). Names should match the values in
#'   \code{species} (or \code{species2} when \code{strains = TRUE}). If
#'   \code{NULL}, a palette is generated by \code{\link{get_palette}}.
#'
#' @details
#' For each feature table, the function:
#' \itemize{
#'   \item Optionally renames strain-level rows using
#'         \code{\link{strain_name2strain_number}} when \code{strains = TRUE}.
#'   \item Filters out features with all-zero counts using
#'         \code{\link{filter_features_by_col_counts}}.
#'   \item Drops samples (columns) that contain only zeros.
#'   \item Converts row names to a \code{species} column and reshapes the table
#'         to long format using \pkg{tidyr::gather}.
#'   \item Adds an \code{experiment} column from \code{experiments_names}.
#' }
#'
#' When \code{strains = TRUE}, additional columns \code{strain} and
#' \code{species2} are created:
#' \itemize{
#'   \item \code{strain}: a label like \code{"Strain 1"}, \code{"Strain 2"}, etc.
#'   \item \code{species2}: the species name without the trailing numeric ID.
#' }
#'
#' The final combined long data frame is filtered to remove zero abundances,
#' converted to factors with \code{experiment} levels matching
#' \code{experiments_names}, and plotted using \pkg{ggplot2}:
#' \itemize{
#'   \item If \code{shared_samples = FALSE}, x-axis = \code{sample}, faceted
#'         by \code{experiment}.
#'   \item If \code{shared_samples = TRUE}, x-axis = \code{experiment}, faceted
#'         by \code{sample}.
#'   \item For \code{strains = TRUE}, \pkg{ggpattern} is used to distinguish
#'         strains by patterns while species are distinguished by fill colors.
#'   \item For \code{strains = FALSE}, a standard stacked barplot is created
#'         with \code{species} as fill.
#' }
#'
#' Relative abundances are shown as fractions of 1 via
#' \code{position = "fill"}.
#'
#' @return A \code{ggplot} object representing the faceted grid of barplots.
#'
#' @examples
#' \dontrun{
#' p <- barplots_grid(
#'   feature_tables   = list(exp1_abund, exp2_abund),
#'   experiments_names = c("Experiment 1", "Experiment 2"),
#'   shared_samples   = FALSE,
#'   strains          = FALSE,
#'   plot_title       = "Relative abundance across experiments"
#' )
#' p
#' }
#'
#' @export
barplots_grid <- function(feature_tables, experiments_names, shared_samples = FALSE, strains = FALSE, plot_title = "",
                          plot_title_size = 14, x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                          y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 0,
                          legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3, legend_key_size = 1,
                          colour_palette = NULL){
  # Creates a grid of Barplots

  ### Step 1. Clean, join and gather the otu tables.
  sample_names = c() # to keep track of the sample names
  for (table in seq(from = 1, to = length(feature_tables), by=1)) { # iterate over all the feature tables
    # copy current feature table to avoid modifying the original table.
    feature_table <- feature_tables[[table]]

    #print(head(feature_table)) # check the working feature table

    if (isTRUE(strains)) {
      # Convert table with strain names to a strain-number table
      feature_table <- strain_name2strain_number(feature_table)
    }

    # Remove rows with Zero counts
    feature_table <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1)

    #print(head(feature_table))

    # save names of species
    species_names <- row.names(feature_table)

    # Remove columns (samples) with zero count
    if (ncol(feature_table) > 1) {
      feature_table <- feature_table[, colSums(feature_table != 0) > 0]
    }

    sample_names <- c(sample_names, colnames(feature_table))

    #print(head(feature_table2))

    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table["species"] <- species_names
    #print(feature_table2$species)

    # Use dplyr gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table, 1:(ncol(feature_table) - 1) , key = "sample", value = "abundance")

    #print(experiments_names[table]) # check experiment name that corresponds to working feature table.

    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector

    #print(head(feature_table_g))

    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data for all tables to make a barplot.
    if (table == 1) {
      plot_df <- feature_table_g
    }else{
      plot_df <- rbind(plot_df, feature_table_g)
    }
  }
  print(sample_names) # check sample_names
  print(head(plot_df)) # check gathered table

  ### Step 2. Convert Strain data to a graphing-compatible format.
  # Add strain data column to long dataframe
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      mutate(
        strain = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)  # Remove strain number from species name
      )
  }

  print(head(plot_df))

  ### Step 3. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    filter(!is.na(abundance) & abundance != 0)

  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      filter(!is.na(strain) & strain != 0)
  }

  plot_df_filtered$experiment <- factor(plot_df_filtered$experiment, levels = experiments_names)

  ### Step 4. Plotting
  # get color palette
  if (is.null(colour_palette)) {
    colour_palette <- get_palette(nColors = length(unique(plot_df$species)))
  }

  print(plot_df_filtered) # check final table prevouos to plotting

  # Create base plot.
  if (shared_samples) {
    p1 <- ggplot(data = plot_df_filtered, aes(x = experiment, y=abundance)) +
      facet_grid(~sample)
  } else{
    p1 <- ggplot(data = plot_df_filtered, aes(x = sample, y=abundance)) +
      facet_grid(~experiment, scales = "free", space = "free")
  }

  # Add elements based on graph type.
  if (isTRUE(strains)) {
    print("strains processing")
    p1 <- p1 + ggpattern::geom_bar_pattern(aes(fill = species2, pattern = strain, pattern_density = strain),
                                           position = "fill",
                                           stat="identity",
                                           show.legend = TRUE,
                                           pattern_color = "white",
                                           pattern_fill = "white",
                                           pattern_angle = 45,
                                           pattern_spacing = 0.025) +
      ggpattern::scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(pattern = guide_legend(override.aes = list(fill = "black")),
             fill = guide_legend(override.aes = list(pattern = "none")))
  } else{
    print("no strains")
    p1 <- p1 + geom_bar(aes(fill = species),
                        position = position_fill(),
                        stat = "identity")
  }

  if (!is.null(colour_palette)) {
    p1 <- p1 + ggplot2::scale_fill_manual(values=colour_palette)
  } else{
    print("Colours vec is null, using standard color palette.")
  }

  p1 <- p1 +
    theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                   axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                   axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, vjust = 0.5, hjust=1, size = x_axis_text_size),
                   axis.title.y = ggplot2::element_text(size=y_axis_title_size, angle = 90),
                   axis.text.y = ggplot2::element_text(size = x_axis_text_size, angle = y_axis_text_angle),
                   legend.title=ggplot2::element_text(size=legend_title_size),
                   legend.text=ggplot2::element_text(size=legend_text_size),
                   legend.position=legend_pos, legend.key.size = unit(legend_key_size, "cm")) +
    guides(fill = guide_legend(ncol = legend_cols))

  # Show plot
  p1

  return(p1)
}


#' Create a Relative Abundance Barplot from a Feature Table
#'
#' Generates a stacked relative abundance barplot from a single feature table
#' (e.g., OTU/ASV or species-by-sample table), with several options for sorting
#' samples and optionally encoding strain-level information using patterns.
#'
#' @param feature_table A numeric matrix or data frame where rows represent
#'   features (e.g., species or strains) and columns represent samples.
#'   Row names are assumed to be feature identifiers.
#' @param sort_type Character string specifying how to order samples on the
#'   x-axis. One of:
#'   \itemize{
#'     \item \code{"none"} (default): keep the original column order.
#'     \item \code{"feature_value"}: order samples by the relative contribution
#'           of a given feature (see \code{feature_to_sort}).
#'     \item \code{"similarity"}: order samples according to hierarchical
#'           clustering via \code{\link{order_samples_by_clustering}}.
#'   }
#' @param feature_to_sort Character string giving the feature (row name) used
#'   for ordering samples when \code{sort_type = "feature_value"}. Ignored
#'   otherwise.
#' @param strains Logical. If \code{TRUE}, treats row names as strain-level
#'   identifiers of the form \code{"Genus species strainID"} and uses
#'   \code{\link{strain_name2strain_number}} plus pattern aesthetics
#'   (\pkg{ggpattern}) to represent strains; species are encoded as fill colors.
#'   Default is \code{FALSE}.
#' @param plot_title Character. Overall plot title. (Note: currently only the
#'   theme title size is controlled; you can extend to use this string explicitly
#'   if desired.) Default is \code{""}.
#' @param plot_title_size Numeric. Text size for the plot title. Default is
#'   \code{14}.
#' @param x_axis_text_size Numeric. Text size for x-axis tick labels. Default
#'   is \code{12}.
#' @param x_axis_title_size Numeric. Text size for the x-axis title. Default
#'   is \code{12}.
#' @param x_axis_text_angle Numeric. Rotation angle (degrees) for x-axis tick
#'   labels. Default is \code{0}.
#' @param y_axis_title_size Numeric. Text size for the y-axis title. Default
#'   is \code{12}.
#' @param y_axis_text_size Numeric. Text size for y-axis tick labels. Default
#*   is \code{12}.
#' @param y_axis_text_angle Numeric. Rotation angle (degrees) for y-axis tick
#'   labels. Default is \code{90}.
#' @param legend_pos Character string specifying the legend position (e.g.,
#'   \code{"right"}, \code{"bottom"}, \code{"none"}). Default is \code{"right"}.
#' @param legend_title_size Numeric. Text size for the legend title. Default is
#'   \code{12}.
#' @param legend_text_size Numeric. Text size for legend text. Default is
#'   \code{12}.
#' @param legend_cols Integer. Number of columns in the legend for the
#'   \code{fill} aesthetic. Default is \code{3}.
#' @param x_vjust Numeric. Vertical justification for x-axis text. Default is
#'   \code{0.5}.
#' @param x_hjust Numeric. Horizontal justification for x-axis text. Default
#'   is \code{1}.
#' @param transform_table Logical. If \code{TRUE} and \code{sort_type =
#'   "similarity"}, the feature table is transformed using
#'   \code{\link{transform_feature_table}} with method \code{"min_max"} prior
#'   to clustering. Default is \code{TRUE}.
#' @param colour_palette Optional named character vector of colors to use for
#'   features (or species groups). If \code{NULL}, a palette is generated with
#'   \code{\link{get_palette}}. When \code{strains = TRUE}, the palette is
#'   applied to \code{species2}.
#' @param replace_c Logical. Passed to \code{\link{get_palette}} as
#'   \code{replace_cols}, controlling whether colors may be sampled with
#'   replacement when generating a palette. Default is \code{FALSE}.
#'
#' @details
#' The function proceeds in several steps:
#' \itemize{
#'   \item Filters out rows (features) with zero counts using
#'         \code{\link{filter_features_by_col_counts}}.
#'   \item Removes samples (columns) that are all zeros.
#'   \item Optionally converts strain-level row names into a standardized
#'         species + numeric strain ID format with
#'         \code{\link{strain_name2strain_number}}.
#'   \item Depending on \code{sort_type}, orders samples by:
#'         \itemize{
#'           \item original column order (\code{"none"}),
#'           \item decreasing relative contribution of a specific feature
#'                 (\code{"feature_value"}),
#'           \item similarity-based clustering
#'                 (\code{"similarity"}, via \code{\link{order_samples_by_clustering}}).
#'         }
#'   \item Reshapes the data into long format with columns \code{species},
#'         \code{sample}, and \code{abundance}. When \code{strains = TRUE},
#'         also creates \code{strain} and \code{species2} columns.
#'   \item Filters out zero and \code{NA} abundances, and, if applicable,
#'         invalid strain entries.
#'   \item Constructs a stacked barplot with relative abundances
#'         (\code{position = "fill"}), using either:
#'         \itemize{
#'           \item \pkg{ggpattern} to encode strains via patterns and species
#'                 via fill colors, or
#'           \item standard stacked bars colored by species.
#'         }
#' }
#'
#' @return A \code{ggplot} object representing the relative abundance barplot.
#'
#' @examples
#' \dontrun{
#' # Basic usage (no sorting)
#' p <- barplot_from_feature_table(feature_table = abundance_mat)
#' p
#'
#' # Sort samples by the relative contribution of a given species
#' p <- barplot_from_feature_table(
#'   feature_table   = abundance_mat,
#'   sort_type       = "feature_value",
#'   feature_to_sort = "Corynebacterium propinquum"
#' )
#'
#' # Sort samples by similarity (clustering-based order)
#' p <- barplot_from_feature_table(
#'   feature_table = abundance_mat,
#'   sort_type     = "similarity",
#'   transform_table = TRUE
#' )
#'
#' # Strain-level visualization with patterns (requires ggpattern)
#' p <- barplot_from_feature_table(
#'   feature_table = strain_level_mat,
#'   strains       = TRUE
#' )
#' }
#'
#' @export
barplot_from_feature_table <- function(feature_table, sort_type = "none", feature_to_sort = NULL, strains = FALSE,
                                       plot_title = "", plot_title_size = 14,
                                       x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                                       y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 90,
                                       legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3,
                                       x_vjust = 0.5, x_hjust = 1, transform_table = TRUE,
                                       colour_palette = NULL, replace_c = FALSE){
  ### Step 1. Clean feature table
  # Remove empty rows (features)
  feature_table2 <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1)

  # Remove columns (samples) with zero count
  if (ncol(feature_table2) > 1) {
    feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
  }

  if (isTRUE(strains)) {
    # Convert table with strain names to a strain-number table
    feature_table2 <- strain_name2strain_number(feature_table2)
  }

  # Save species names from row names
  species <- row.names(feature_table2)
  print(head(feature_table2))

  ### Step 2. If sorting, determine sample order.
  if (sort_type == "feature_value" && !is.null(feature_to_sort)) {
    print("Sort samples by feature_value")
    # Make "species" column with the rownames
    df1 <- feature_table2 %>% rownames_to_column(var = "species")

    total_abundance <- colSums(df1[, -1])

    # Filter the row of the species of interest and calculate its proportion with respect to total abundance
    df_proportion <- df1 %>%
      dplyr::filter(species == feature_to_sort) %>%
      dplyr::select(-species)

    # calculate species of interest proportion
    df_proportion <- df_proportion[1, ] / total_abundance

    # Get sample names sorted by the species of interest proportion
    ordered_samples <- df_proportion %>%
      unlist() %>%
      sort(decreasing = TRUE) %>%
      names()

  } else if (sort_type == "similarity") {
    print("Sort samples by similarity")

    # transform table
    if (transform_table) {
      df1 <- transform_feature_table(feature_table = feature_table2, transform_method = "min_max")
    } else {
      df1 <- feature_table2
    }

    # Get the order of samples based on clustering
    ordered_samples <- order_samples_by_clustering(df1)

    df1 <- df1 %>% rownames_to_column(var = "species")

  } else if (sort_type == "none") {
    print("No sorting chosen")
    df1 <- feature_table2
    ordered_samples <- colnames(feature_table2)
    # Generate a column with the names of ASVs/OTUs using rownames.
    df1["species"] <- species
  } else {
    print("No valid sorting option chosen")
    return()
  }

  ### Step 3. Process features table to plotting table.
  plot_df <- df1 %>%
    tidyr::pivot_longer(-species, names_to = "sample", values_to = "abundance")

  # If strain processing has to be done.
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      dplyr::mutate(
        strain   = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)                   # Remove strain number from species name
      )
  }

  ### Step 4. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    dplyr::filter(!is.na(abundance) & abundance != 0)

  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      dplyr::filter(!is.na(strain) & strain != 0)
  }

  # Factor the "sample" variable so the order of samples is as in "ordered_samples"
  plot_df_filtered$sample <- factor(plot_df_filtered$sample, levels = ordered_samples)

  print(head(plot_df_filtered))

  ### Step 5. Create plot.
  if (is.null(colour_palette)) { # get colour palette
    print("Colour pallette generated")
    nfeatures <- length(unique(plot_df_filtered$species))
    colour_palette <- get_palette(nColors = nfeatures, replace_cols = replace_c)
    print(colour_palette)
  }

  # Create base plot.
  ft_barplot <- ggplot2::ggplot(plot_df_filtered, ggplot2::aes(x = sample, y = abundance, fill = species))

  if (isTRUE(strains)) {
    print("strains processing")
    ft_barplot <- ft_barplot +
      ggpattern::geom_bar_pattern(
        ggplot2::aes(fill = species2, pattern = strain, pattern_density = strain),
        position        = "fill",
        stat            = "identity",
        show.legend     = TRUE,
        pattern_color   = "white",
        pattern_fill    = "white",
        pattern_angle   = 45,
        pattern_spacing = 0.025
      ) +
      ggpattern::scale_pattern_manual(
        values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")
      ) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(
        pattern = guide_legend(override.aes = list(fill = "black")),
        fill    = guide_legend(override.aes = list(pattern = "none"))
      )
  } else {
    print("no strains")
    ft_barplot <- ft_barplot +
      ggplot2::geom_bar(
        ggplot2::aes(fill = species),
        position = ggplot2::position_fill(),
        stat     = "identity"
      )
  }

  # add theme options
  ft_barplot <- ft_barplot +
    ggplot2::theme_void() +
    ggplot2::scale_fill_manual(values = colour_palette) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = plot_title_size, face = "bold"),
      axis.title.x = ggplot2::element_text(size = x_axis_title_size),
      axis.text.x  = ggplot2::element_text(
        angle = x_axis_text_angle,
        vjust = x_vjust,
        hjust = x_hjust,
        size  = x_axis_text_size
      ),
      axis.title.y = ggplot2::element_text(
        size   = y_axis_title_size,
        angle  = y_axis_text_angle,
        margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)
      ),
      axis.text.y  = ggplot2::element_text(size = y_axis_text_size),
      legend.position = legend_pos,
      legend.title    = ggplot2::element_text(size = legend_title_size),
      legend.text     = ggplot2::element_text(size = legend_text_size)
    ) +
    guides(fill = guide_legend(ncol = legend_cols))

  ft_barplot
}


# ----- PCoA -----
#' Align Metabolomics Feature Table with Metadata Attributes
#'
#' Aligns the columns (samples) of a metabolomics feature table to a metadata
#' table and ensures that key attribute columns needed for downstream analyses
#' are present and correctly typed.
#'
#' The function identifies a sample ID column in the metadata (or derives it
#' from row names), intersects sample IDs between the metabolomics table and
#' metadata, reorders both objects to the same sample order, and optionally
#' derives missing attributes from the sample ID strings.
#'
#' @param metab_df A numeric matrix or data frame where columns are samples and
#'   rows are features (e.g., metabolites). Column names must contain sample IDs.
#' @param metadata_df A data frame containing sample metadata. Sample IDs must
#'   be available either in a column (see \code{sample_col}) or in row names.
#' @param sample_col Optional character string. Name of the column in
#'   \code{metadata_df} containing sample IDs. If \code{NULL} (default), the
#'   function attempts to use a \code{"Sample"} column if present, otherwise
#'   uses row names as sample IDs and creates a \code{"Sample"} column.
#'
#' @details
#' The function expects or constructs the following metadata columns:
#' \itemize{
#'   \item \code{ATTRIBUTE_Cluster}: cluster assignment (factor).
#'   \item \code{ATTRIBUTE_Time}: timepoint information. If missing, it is
#'         derived from the sample ID using the pattern \code{"_T<digits>_"}
#'         (e.g., \code{"SC01_T3_rep1"} gives timepoint \code{3}).
#'   \item \code{ATTRIBUTE_SynCom}: SynCom identifier. If missing, it is derived
#'         as the substring before the first underscore in the sample ID
#'         (e.g., \code{"SC01_T3_rep1"} gives \code{"SC01"}).
#' }
#'
#' Only overlapping sample IDs between \code{colnames(metab_df)} and the
#' metadata sample ID column are retained. The metabolomics matrix is converted
#' to a numeric matrix and returned alongside the aligned metadata.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{X}{A numeric matrix containing \code{metab_df} restricted to shared
#'     samples and ordered to match the metadata. Columns are samples.}
#'   \item{meta}{A data frame containing \code{metadata_df} restricted to shared
#'     samples and ordered to match \code{X}. Row names are set to sample IDs.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- align_samples_attr(
#'   metab_df    = metab_mat,
#'   metadata_df = sample_metadata,
#'   sample_col  = "Sample"
#' )
#'
#' X  <- res$X
#' md <- res$meta
#' }
#'
#' @export
align_samples_attr <- function(metab_df, metadata_df, sample_col = NULL) {
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df

  # Figure out where sample IDs live in metadata
  if (is.null(sample_col)) {
    if ("Sample" %in% colnames(md)) {
      sample_col <- "Sample"
    } else if (!is.null(rownames(md))) {
      md$Sample <- rownames(md)
      sample_col <- "Sample"
    } else stop("metadata_df must have a 'Sample' column or rownames = sample IDs.")
  }

  # Ensure required columns exist; fill from IDs if missing
  req <- c("ATTRIBUTE_Cluster", "ATTRIBUTE_Time", "ATTRIBUTE_SynCom")
  missing <- setdiff(req, colnames(md))
  if (length(missing) > 0) {
    if ("ATTRIBUTE_Cluster" %in% missing) md$ATTRIBUTE_Cluster <- NA
    if ("ATTRIBUTE_Time" %in% missing) {
      md$ATTRIBUTE_Time <- suppressWarnings(
        as.integer(sub(".*_T(\\d+)_.*", "\\1", md[[sample_col]]))
      )
    }
    if ("ATTRIBUTE_SynCom" %in% missing) {
      md$ATTRIBUTE_SynCom <- sub("_.*", "", md[[sample_col]])  # "SC##"
    }
  }

  # Intersect & align
  common <- intersect(colnames(metab_df), md[[sample_col]])
  if (length(common) == 0) stop("No overlapping sample IDs between metabolome and metadata.")

  X  <- as.matrix(metab_df[, common, drop = FALSE])
  md <- md[match(common, md[[sample_col]]), , drop = FALSE]
  rownames(md) <- md[[sample_col]]

  # Types
  md$ATTRIBUTE_Cluster <- factor(md$ATTRIBUTE_Cluster)
  md$ATTRIBUTE_Time    <- factor(md$ATTRIBUTE_Time)
  md$ATTRIBUTE_SynCom  <- as.character(md$ATTRIBUTE_SynCom)

  print("Data alligned succesfully")

  list(X = X, meta = md)
}


#' Flexible PCoA workflow with plotting and optional PERMANOVA
#'
#' Computes a Principal Coordinates Analysis (PCoA) from a feature table
#' (e.g., metabolomics) and sample metadata, returning ordination results,
#' variance explained, plotting-ready scores, a ggplot object, the distance
#' object/matrix, and an optional single-factor PERMANOVA.
#'
#' This function:
#' \itemize{
#'   \item aligns samples between \code{metab_df} columns and \code{metadata_df}
#'   \item optionally preprocesses data (none / sqrt / hellinger)
#'   \item computes Bray-Curtis or Euclidean distances
#'   \item runs \code{cmdscale()} PCoA (with \code{add = TRUE})
#'   \item makes a \code{ggplot2} ordination plot with optional ellipses and labels
#'   \item optionally runs \code{vegan::adonis2()} PERMANOVA for one variable
#' }
#'
#' @param metab_df Feature-by-sample table. Columns must be sample IDs. Values
#'   are numeric (e.g., metabolite intensities/abundances).
#' @param metadata_df Sample metadata. Sample IDs can be provided via rownames
#'   (default) or via \code{sample_id_col}.
#' @param sample_id_col Optional. Column name in \code{metadata_df} that contains
#'   sample IDs. If \code{NULL}, rownames(\code{metadata_df}) are used as IDs.
#' @param color_var Metadata column name used to color points (treated as factor).
#' @param shape_var Optional. Metadata column name used to shape points (factor).
#' @param ellipse_var Optional. Metadata column name used to define ellipse groups
#'   (factor). Ellipses are drawn only if \code{ellipse = TRUE}.
#' @param distance Distance metric. One of \code{"bray"} or \code{"euclidean"}.
#' @param preprocess Preprocessing applied before distance calculation. One of
#'   \code{"none"}, \code{"hellinger"}, or \code{"sqrt"}.
#'   \itemize{
#'     \item \code{"hellinger"} uses \code{vegan::decostand(method = "hellinger")}
#'     \item \code{"sqrt"} applies \code{sqrt()}
#'   }
#' @param k_axes Number of PCoA axes to compute (minimum 2 are computed internally).
#' @param ellipse Logical. If \code{TRUE} and \code{ellipse_var} is provided,
#'   draw 95% normal ellipses via \code{ggplot2::stat_ellipse()}.
#' @param min_n_for_ellipse Minimum group size required (per \code{ellipse_var}
#'   level) to draw ellipses. Groups with fewer samples are skipped.
#' @param label_points Logical. If \code{TRUE}, label points with sample IDs using
#'   \code{ggrepel::geom_text_repel()}.
#' @param points_palette Optional named character vector mapping
#'   \code{levels(metadata_df[[color_var]])} to colors, passed to
#'   \code{scale_color_manual()}.
#' @param ellipse_palette Optional named character vector mapping
#'   \code{levels(metadata_df[[ellipse_var]])} to colors when ellipses use a
#'   separate color scale.
#' @param color_var_leg_columns Integer. Number of legend columns for the color
#'   legend (points).
#' @param permanova_var Optional. Metadata column to test in PERMANOVA (single-term).
#'   Internally mapped to \code{.__var} and tested as \code{dd ~ .__var}.
#' @param strata_var Optional. Metadata column used as a blocking factor (strata)
#'   in PERMANOVA (e.g., SynCom).
#' @param permutations Integer. Number of permutations for PERMANOVA.
#'
#' @details
#' **Sample alignment**
#' Samples are intersected between \code{colnames(metab_df)} and metadata sample
#' IDs. The feature table is subset to overlapping samples and metadata is
#' re-ordered to match the feature table column order.
#'
#' **Distance and preprocessing constraints**
#' \itemize{
#'   \item Bray-Curtis requires non-negative values; the function errors if any
#'     negative entries are present.
#'   \item Hellinger and sqrt preprocessing require non-negative values.
#' }
#'
#' **Ellipses**
#' Ellipses are drawn only for groups with at least \code{min_n_for_ellipse}
#' samples. If \code{ellipse_var != color_var}, a second color scale is created
#' using \code{ggnewscale::new_scale_color()} (loaded on demand).
#'
#' **PERMANOVA**
#' If \code{permanova_var} is provided, rows with missing values in the tested
#' variable (and in \code{strata_var}, if provided) are dropped before running
#' \code{vegan::adonis2()}. The distance object is subset accordingly.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{pcoa}: \code{cmdscale()} result (includes eigenvalues and points)
#'   \item \code{explained}: percent variance explained for positive eigenvalues
#'   \item \code{scores}: data.frame with PCoA coordinates + metadata + Sample
#'   \item \code{plot}: ggplot object of PCoA ordination
#'   \item \code{dist}: \code{vegan::vegdist} distance object
#'   \item \code{dist_mat}: distance matrix as a base \code{matrix}
#'   \item \code{permanova}: \code{adonis2} result or \code{NULL}
#'   \item \code{settings}: list of key settings used for reproducibility
#' }
#'
#' @examples
#' \dontrun{
#' res <- pcoa_flex(
#'   metab_df = feat_table,
#'   metadata_df = md,
#'   color_var = "Cluster",
#'   shape_var = "Time",
#'   ellipse_var = "Cluster",
#'   distance = "bray",
#'   preprocess = "hellinger",
#'   permanova_var = "Cluster",
#'   strata_var = "SynCom"
#' )
#' res$plot
#' res$permanova
#' }
#'
#' @export
pcoa_flex <- function(
    metab_df, metadata_df,
    sample_id_col = NULL,              # "Sample" if IDs live in a column; NULL -> use rownames(metadata_df)
    color_var,                          # metadata column for point color (factor)
    shape_var = NULL,                   # metadata column for point shape (factor) or NULL (Optional)
    ellipse_var = NULL,                 # metadata column for ellipses (factor) or NULL (Optional)
    distance = c("bray","euclidean"),
    preprocess = c("none","hellinger","sqrt"),
    k_axes = 2,
    ellipse = TRUE,
    min_n_for_ellipse = 3,
    label_points = FALSE,
    points_palette  = NULL,             # named vector for levels(color_var) (optional)
    ellipse_palette = NULL,             # named vector for levels(ellipse_var) (optional)
    color_var_leg_columns = 1,
    # PERMANOVA options (single variable test)
    permanova_var = NULL,               # metadata column to test (factor or numeric)
    strata_var = NULL,                  # optional blocking factor (e.g., SynCom)
    permutations = 999
) {
  distance   <- match.arg(distance)
  preprocess <- match.arg(preprocess)

  # ---- Load required packages (quietly) ----
  # Loads ggplot2 + vegan if not already attached. (ggnewscale and ggrepel are
  # loaded on demand later depending on options.)
  pkgs <- c("ggplot2","vegan")
  to_load <- pkgs[!pkgs %in% (.packages())]
  if (length(to_load)) suppressPackageStartupMessages(sapply(to_load, require, character.only = TRUE))

  # ---------- Align samples ----------
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df

  # Find sample IDs in metadata (column or rownames)
  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(md)) stop("sample_id_col not found in metadata.")
    ids <- as.character(md[[sample_id_col]])
  } else {
    if (is.null(rownames(md))) stop("metadata_df must have rownames or provide sample_id_col.")
    ids <- rownames(md)
  }

  # Restrict to common sample IDs and enforce consistent ordering
  common <- intersect(colnames(metab_df), ids)
  if (!length(common)) stop("No overlapping sample IDs between matrix columns and metadata IDs.")

  # Subset the feature table to shared samples and convert to matrix
  X  <- as.matrix(metab_df[, common, drop = FALSE]) # features x samples

  # Reorder metadata rows to match X column order
  md <- if (!is.null(sample_id_col)) md[match(common, ids), , drop = FALSE]
  else                              md[match(common, rownames(md)), , drop = FALSE]
  rownames(md) <- common

  # Validate and coerce aesthetic variables to factors (dropping unused levels)
  if (!color_var %in% colnames(md))   stop("color_var not found in metadata.")
  if (!is.null(shape_var)   && !shape_var   %in% colnames(md)) stop("shape_var not found in metadata.")
  if (!is.null(ellipse_var) && !ellipse_var %in% colnames(md)) stop("ellipse_var not found in metadata.")

  md[[color_var]] <- droplevels(factor(md[[color_var]]))
  if (!is.null(shape_var))   md[[shape_var]]   <- droplevels(factor(md[[shape_var]]))
  if (!is.null(ellipse_var)) md[[ellipse_var]] <- droplevels(factor(md[[ellipse_var]]))

  # ---------- Build samples x features and preprocess ----------
  # cmdscale/vegdist typically operate on samples as rows
  S <- t(X)  # samples x metabolites

  # Distance calculation with safeguards for method requirements
  if (distance == "bray") {
    # Bray-Curtis requires non-negative values
    if (any(S < 0, na.rm = TRUE))
      stop("Bray-Curtis requires non-negative data. Use a non-negative matrix or switch distance='euclidean'.")

    # Optional preprocessing for compositional/abundance-style data
    if (preprocess == "hellinger")      S <- vegan::decostand(S, method = "hellinger")
    else if (preprocess == "sqrt")      S <- sqrt(S)

    dd <- vegan::vegdist(S, method = "bray")
  } else { # euclidean
    # Allow preprocessing, but enforce non-negativity if the chosen transform needs it
    if (preprocess %in% c("hellinger","sqrt")) {
      if (any(S < 0, na.rm = TRUE)) stop("Selected preprocessing requires non-negative data.")
      if (preprocess == "hellinger") S <- vegan::decostand(S, method = "hellinger")
      if (preprocess == "sqrt")      S <- sqrt(S)
    }
    dd <- vegan::vegdist(S, method = "euclidean")
  }

  # ---------- PCoA ----------
  # Compute at least 2 axes (required for plotting) and up to k_axes axes returned
  pc <- cmdscale(dd, eig = TRUE, k = max(2, k_axes), add = TRUE)

  # Assemble coordinates + metadata into a plotting-ready data.frame
  scores <- as.data.frame(pc$points)
  colnames(scores) <- paste0("PCo", seq_len(ncol(scores)))
  scores$Sample <- rownames(scores)
  scores <- cbind(scores, md[rownames(scores), , drop = FALSE])

  # Percent variance explained (computed using positive eigenvalues only)
  eig <- pc$eig
  pos <- eig[eig > 0]
  explained <- 100 * pos / sum(pos)

  # ---------- Plot (tidy evaluation; no aes_string) ----------
  # Build dynamic point aesthetics (color always; shape optional)
  pt_map <- if (is.null(shape_var)) {
    ggplot2::aes(color = .data[[color_var]])
  } else {
    ggplot2::aes(color = .data[[color_var]], shape = .data[[shape_var]])
  }

  # Base ordination scatter plot
  p <- ggplot2::ggplot(scores, ggplot2::aes(PCo1, PCo2)) +
    ggplot2::geom_point(mapping = pt_map, size = 2, alpha = 0.9) +
    ggplot2::labs(
      x = paste0("PCo1 (", round(explained[1], 1), "%)"),
      y = paste0("PCo2 (", round(explained[2], 1), "%)"),
      color = color_var,
      shape = if (!is.null(shape_var)) shape_var else NULL
    ) +
    ggplot2::theme_bw()

  # Optional manual palette for point colors, otherwise use discrete scale
  if (!is.null(points_palette)) {
    p <- p + ggplot2::scale_color_manual(
      values = points_palette,
      name = color_var,
      guide = guide_legend(ncol = color_var_leg_columns)
    )
  } else {
    p <- p + ggplot2::scale_color_discrete(
      name = color_var,
      guide = guide_legend(ncol = color_var_leg_columns)
    )
  }

  # ---------- Ellipses (optional) ----------
  # Draw ellipses only if enabled and grouping variable provided, and only for
  # groups with at least min_n_for_ellipse samples.
  if (isTRUE(ellipse) && !is.null(ellipse_var)) {
    grp_counts <- table(scores[[ellipse_var]])
    ok_groups  <- names(grp_counts)[grp_counts >= min_n_for_ellipse]

    if (length(ok_groups)) {
      dat_ell <- scores[scores[[ellipse_var]] %in% ok_groups, , drop = FALSE]

      if (ellipse_var == color_var) {
        # Reuse the same color scale as points
        p <- p + ggplot2::stat_ellipse(
          data = dat_ell,
          ggplot2::aes(
            x = PCo1, y = PCo2,
            group = .data[[ellipse_var]],
            color = .data[[ellipse_var]]
          ),
          inherit.aes = FALSE,
          level = 0.95, type = "norm",
          linewidth = 0.6, alpha = 0.9
        )
      } else {
        # Use a second color scale for ellipses if ellipse grouping differs from point color
        #if (!"ggnewscale" %in% (.packages())) suppressPackageStartupMessages(require(ggnewscale))

        p <- p +
          ggnewscale::new_scale_color() +
          ggplot2::stat_ellipse(
            data = dat_ell,
            ggplot2::aes(
              x = PCo1, y = PCo2,
              group = .data[[ellipse_var]],
              color = .data[[ellipse_var]]
            ),
            inherit.aes = FALSE,
            level = 0.95, type = "norm",
            linewidth = 0.7, alpha = 0.9
          ) +
          {
            if (!is.null(ellipse_palette))
              ggplot2::scale_color_manual(
                values = ellipse_palette,
                name = paste0(ellipse_var, " (ellipse)")
              )
            else
              ggplot2::scale_color_discrete(
                name = paste0(ellipse_var, " (ellipse)")
              )
          }
      }
    } else {
      message(
        "Skipping ellipses: fewer than ", min_n_for_ellipse,
        " samples in all groups of '", ellipse_var, "'."
      )
    }
  }

  # ---------- Labels on points (optional) ----------
  if (label_points) {
    #if (!"ggrepel" %in% (.packages())) suppressPackageStartupMessages(require(ggrepel))
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = Sample),
      size = 2,
      max.overlaps = 60
    )
  }

  # ---------- PERMANOVA (single variable; optional) ----------
  permanova <- NULL
  if (!is.null(permanova_var)) {
    if (!permanova_var %in% colnames(md)) stop("permanova_var not found in metadata.")

    # Keep only complete cases for the tested variable (and strata if used)
    keep <- complete.cases(md[, permanova_var, drop = FALSE]) &
      (if (!is.null(strata_var)) complete.cases(md[, strata_var, drop = FALSE]) else TRUE)

    if (!all(keep)) message("PERMANOVA: dropping ", sum(!keep), " sample(s) with NA in selected variables.")
    md_perm <- md[keep, , drop = FALSE]

    # Subset the distance object to the kept samples
    labs <- rownames(md_perm)
    dd_mat <- as.matrix(dd)
    dd_perm <- stats::as.dist(dd_mat[labs, labs])

    # Build one-term model dd_perm ~ .__var
    md_perm$.__var <- md_perm[[permanova_var]]
    form <- stats::as.formula("dd_perm ~ .__var")

    permanova <- vegan::adonis2(
      formula = form,
      data = md_perm,
      permutations = permutations,
      strata = if (!is.null(strata_var)) md_perm[[strata_var]] else NULL
    )
  }

  # Return results + useful objects for downstream use
  list(
    pcoa      = pc,
    explained = explained,
    scores    = scores,
    plot      = p,
    dist      = dd,
    dist_mat  = as.matrix(dd),
    permanova = permanova,
    settings  = list(
      distance = distance, preprocess = preprocess,
      color_var = color_var, shape_var = shape_var,
      ellipse_var = ellipse_var, permanova_var = permanova_var,
      strata_var = strata_var, permutations = permutations
    )
  )
}


# ----- Targeted metabolomics analyses -----
#' Extract replicate-aware sample information from a data frame
#'
#' Identifies columns that look like replicate measurements (e.g., \code{Sample_1},
#' \code{Sample_2}, ...) using a regular expression, then derives the base sample
#' name (prefix) for each replicate column and the set of unique base samples.
#'
#' This is useful when your wide table stores technical/biological replicates as
#' separate columns that share a common prefix and differ only by a trailing
#' suffix after the final underscore.
#'
#' @param df A data frame with replicate columns in its column names.
#' @param replicate_regex Character string. Regular expression used to identify
#'   replicate columns. The default (\code{"^[^_]+_\\\\d+$"}) matches names with
#'   a single underscore and a trailing integer (e.g., \code{"ABC_1"}). If your
#'   naming includes additional underscores (e.g., \code{"Subject_A_1"}), provide
#'   a different pattern.
#'
#' @details
#' Steps performed:
#' \itemize{
#'   \item Find column names matching \code{replicate_regex}.
#'   \item Strip the final underscore segment to define the base sample name
#'     (e.g., \code{"Sample_2"} \eqn{\rightarrow} \code{"Sample"}).
#'   \item Compute unique base sample names.
#' }
#'
#' @return A list with:
#' \itemize{
#'   \item \code{sample_cols}: character vector of replicate column names.
#'   \item \code{base_names}: named character vector mapping each replicate
#'     column to its base sample name (names are \code{sample_cols}).
#'   \item \code{unique_samples}: character vector of unique base sample names.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Met = c("M1","M2"),
#'   A_1 = c(1,2),
#'   A_2 = c(3,4),
#'   B_1 = c(5,6)
#' )
#' info <- get_sample_info(df)
#' info$sample_cols
#' info$base_names
#' info$unique_samples
#'
#' # If sample names can contain underscores (e.g., "Subject_A_1"):
#' info2 <- get_sample_info(df, replicate_regex = "^.+_\\\d+$")
#' }
#'
#' @export
get_sample_info <- function(df, replicate_regex = "^[^_]+_\\d+$") {
  # 1) Find replicate columns like PREFIX_1, PREFIX_2, ...
  sample_cols <- grep(replicate_regex, colnames(df), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No replicate columns found. Check your column names and 'replicate_regex'.")
  }

  # 2) Extract base names (prefix before final underscore)
  # E.g., "Sample_1" -> "Sample"
  base_names <- sub("_[^_]+$", "", sample_cols)

  # 3) Unique sample prefixes
  unique_samples <- unique(base_names)

  # 4) Return everything useful
  list(
    sample_cols    = sample_cols,
    base_names     = setNames(base_names, sample_cols), # map replicate col -> base name
    unique_samples = unique_samples
  )
}


#' Build raw and replicate-mean matrices from a wide data frame
#'
#' Given a data frame containing metabolite rows and replicate sample columns
#' (e.g., \code{Sample_1}, \code{Sample_2}, ...), this function:
#' \itemize{
#'   \item determines metabolite row identifiers (rownames or \code{Metabolite} column)
#'   \item extracts replicate columns into a numeric matrix (\code{mat_raw})
#'   \item computes per-sample means across replicates based on provided prefixes
#'     (\code{mat_mean})
#' }
#'
#' @param df A data frame containing metabolite identifiers and replicate columns.
#' @param sample_cols Character vector of replicate column names to extract from
#'   \code{df} (e.g., those returned by \code{get_sample_info()}).
#' @param base_names Named character vector mapping replicate column names to
#'   base sample names (prefixes). Names should be replicate column names, and
#'   values are their corresponding base sample IDs.
#'
#' @details
#' **Metabolite identifiers**
#' \itemize{
#'   \item If \code{rownames(df)} are present, they are used.
#'   \item Otherwise, if a \code{Metabolite} column exists, it is used.
#'   \item Missing/empty names are replaced by \code{"NA_metabolite"}.
#'   \item Duplicate metabolite names are made unique via \code{make.unique()}.
#' }
#'
#' **Replicate means**
#' Replicates are grouped by \code{base_names} (prefix). For each unique base
#' sample, the function computes \code{rowMeans(..., na.rm = TRUE)} across its
#' replicate columns.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{mat_raw}: numeric matrix of metabolites x replicate-columns.
#'   \item \code{mat_mean}: numeric matrix of metabolites x unique base samples
#'     (replicate means).
#'   \item \code{unique_samples}: character vector of unique base sample IDs.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Metabolite = c("M1","M2"),
#'   A_1 = c(1, 2),
#'   A_2 = c(3, 4),
#'   B_1 = c(5, 6)
#' )
#' si <- get_sample_info(df)
#' mats <- build_mats_from_df(df, si$sample_cols, si$base_names)
#' mats$mat_raw
#' mats$mat_mean
#' }
#'
#' @export
build_mats_from_df <- function(df, sample_cols, base_names) {
  # ---- 0) Determine rownames for metabolites ----
  # Prefer existing rownames(df). Otherwise, fall back to a "Metabolite" column.
  if (is.null(rownames(df))) {
    if ("Metabolite" %in% names(df)) {
      rn <- as.character(df$Metabolite)
    } else {
      stop("No rownames and no 'Metabolite' column found to set rownames.")
    }
  } else {
    rn <- rownames(df)
  }

  # Normalize missing/blank metabolite IDs and enforce uniqueness
  rn[is.na(rn) | rn == ""] <- "NA_metabolite"
  if (anyDuplicated(rn)) rn <- make.unique(rn)

  # ---- 1) Pull replicate columns and coerce to a numeric matrix ----
  if (length(sample_cols) == 0) stop("sample_cols is empty.")
  X <- df[, sample_cols, drop = FALSE]

  # Guard against accidental non-numeric columns (coerce safely)
  if (!all(vapply(X, is.numeric, logical(1)))) {
    X <- data.matrix(X)
  } else {
    X <- as.matrix(X)
  }
  rownames(X) <- rn

  # ---- 2) Compute per-sample (prefix) means across replicates ----
  # base_names should be a named vector: names = replicate cols, values = prefixes
  if (is.null(names(base_names))) {
    # If names are missing, assume base_names corresponds to sample_cols order
    names(base_names) <- sample_cols
  }

  unique_samples <- unique(unname(base_names))

  # For each base sample, average across its replicate columns (NA-safe)
  mat_mean <- sapply(unique_samples, function(smpl) {
    cols <- names(base_names)[base_names == smpl]
    rowMeans(X[, cols, drop = FALSE], na.rm = TRUE)
  })
  mat_mean <- as.matrix(mat_mean)
  rownames(mat_mean) <- rownames(X)

  list(
    mat_raw        = X,
    mat_mean       = mat_mean,
    unique_samples = unique_samples
  )
}



#' Compute log2 fold-changes vs control and significance "stars" (FDR + effect-size gated)
#'
#' Computes (i) a log2 fold-change (LFC) matrix comparing each treatment sample
#' (prefix) against a specified control prefix using replicate means, and (ii) a
#' corresponding matrix of significance markers ("*") based on replicate-level
#' Welch t-tests with Benjamini-Hochberg FDR correction, additionally gated by
#' a minimum absolute LFC threshold.
#'
#' The intent is to match what is visualized in a heatmap:
#' \itemize{
#'   \item LFC values are computed from \code{mat_mean} (means across replicates),
#'     using \code{pseudocount_disp}.
#'   \item Significance is computed from replicate-level data in \code{mat_raw}
#'     after log2-transform with \code{pseudocount_test}, then BH-adjusted, and
#'     only called significant if both FDR and effect-size gates are passed.
#' }
#'
#' @param mat_raw Numeric matrix of metabolites x replicate-columns (raw replicate values).
#'   Row names must identify metabolites.
#' @param mat_mean Numeric matrix of metabolites x sample-prefix columns (replicate means).
#'   Column names are expected to be the unique sample prefixes (including control).
#' @param base_names Named character vector mapping replicate column names to base sample
#'   prefixes (names = replicate columns, values = prefixes). Used to select control and
#'   treatment replicate columns in \code{mat_raw}.
#' @param control_prefix Character. The sample prefix in \code{mat_mean} (and \code{base_names})
#'   that defines the control condition (default \code{"CTRL"}).
#' @param alpha Numeric. FDR cutoff applied to BH-adjusted p-values (default \code{0.05}).
#' @param lfc_gate Numeric. Minimum absolute LFC required to mark significance
#'   (default \code{2}, i.e., 4-fold on the original scale).
#' @param pseudocount_test Numeric. Pseudocount added to replicate intensities in \code{mat_raw}
#'   before log2 transform for t-tests. Helps handle zeros (default \code{1}).
#' @param pseudocount_disp Numeric. Pseudocount added to means in \code{mat_mean} when computing
#'   the heatmap LFC. Kept small to minimally perturb ratios (default \code{1e-8}).
#'
#' @details
#' **LFC definition (heatmap)**
#' For each metabolite \eqn{m} and treatment sample \eqn{s}:
#' \deqn{LFC_{m,s} = \log_2\left(\frac{\bar{x}_{m,s} + \epsilon}{\bar{x}_{m,ctrl} + \epsilon}\right)}
#' where \eqn{\bar{x}} are replicate means from \code{mat_mean} and \eqn{\epsilon = pseudocount_disp}.
#' The control column is dropped from the returned LFC matrix.
#'
#' **Significance testing**
#' For each treatment prefix \code{smpl}, replicate columns are selected via \code{base_names}.
#' Values are transformed using \code{log2(value + pseudocount_test)} and compared against the
#' corresponding control replicates with Welch's t-test (\code{var.equal = FALSE}).
#' BH correction is applied across metabolites for each treatment-vs-control comparison.
#'
#' **Gating**
#' A metabolite is marked with "*" for a given treatment if:
#' \itemize{
#'   \item \code{padj < alpha}, and
#'   \item \code{abs(LFC_heatmap) >= lfc_gate}
#' }
#' where \code{LFC_heatmap} is taken from the same LFC matrix used for plotting.
#'
#' Metabolites with insufficient non-missing replicate values (<2 per group) yield \code{NA}
#' p-values and are never starred.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{lfc}: numeric matrix of metabolites x treatment prefixes (control removed).
#'   \item \code{stars}: character matrix of same dimensions as \code{lfc}, containing "*" for
#'   significant entries and "" otherwise.
#' }
#'
#' @examples
#' \dontrun{
#' # mat_raw: metabolites x replicate columns (e.g., A_1, A_2, CTRL_1, CTRL_2)
#' # mat_mean: metabolites x prefixes (e.g., A, CTRL)
#' res <- compute_lfc_and_stars(
#'   mat_raw = mat_raw,
#'   mat_mean = mat_mean,
#'   base_names = base_names,
#'   control_prefix = "CTRL",
#'   alpha = 0.05,
#'   lfc_gate = 1
#' )
#' res$lfc
#' res$stars
#' }
#'
#' @export
compute_lfc_and_stars <- function(mat_raw,
                                  mat_mean,
                                  base_names,             # named vector: replicate_col -> prefix
                                  control_prefix = "CTRL",
                                  alpha = 0.05,           # FDR cutoff (BH)
                                  lfc_gate = 2,           # |log2 FC| threshold (2 = 4x)
                                  pseudocount_test = 1,   # added to replicate intensities before log2 for t-tests
                                  pseudocount_disp = 1e-8 # added to means before building heatmap LFC
) {
  # ---- Input checks ----
  stopifnot(!is.null(rownames(mat_raw)), !is.null(rownames(mat_mean)))

  unique_samples <- colnames(mat_mean)
  if (is.null(unique_samples)) stop("mat_mean must have column names (sample prefixes).")
  if (!control_prefix %in% unique_samples)
    stop(sprintf("control_prefix '%s' not found in mat_mean columns.", control_prefix))

  # --- 1) Heatmap LFC vs control (means across replicates) ---
  # LFC is computed from the same matrix that would be visualized in the heatmap.
  ctrl_mean <- mat_mean[, control_prefix, drop = TRUE]
  lfc_heatmap_full <- log2((mat_mean + pseudocount_disp) / (ctrl_mean + pseudocount_disp))

  # Drop control column for plotting/stars
  lfc <- lfc_heatmap_full[, setdiff(colnames(lfc_heatmap_full), control_prefix), drop = FALSE]

  # --- 2) Prepare stars matrix (same shape as lfc) ---
  stars <- matrix(
    "",
    nrow = nrow(lfc),
    ncol = ncol(lfc),
    dimnames = dimnames(lfc)
  )

  # --- 3) T-test on replicate-level log2 values ---
  # Use replicate columns defined by base_names, comparing each treatment prefix to control.
  rep_cols  <- names(base_names)
  ctrl_cols <- rep_cols[base_names == control_prefix]
  log_ctrl  <- log2(mat_raw[, ctrl_cols, drop = FALSE] + pseudocount_test)

  for (smpl in colnames(lfc)) {
    trt_cols <- rep_cols[base_names == smpl]
    if (length(trt_cols) == 0) next

    log_trt <- log2(mat_raw[, trt_cols, drop = FALSE] + pseudocount_test)

    # Compute p-values per metabolite (Welch t-test), guarding against insufficient data
    pvals <- vapply(seq_len(nrow(mat_raw)), function(i) {
      x <- log_trt[i, ]
      y <- log_ctrl[i, ]
      if (all(is.na(x)) || all(is.na(y)) ||
          length(na.omit(x)) < 2 || length(na.omit(y)) < 2) return(NA_real_)
      tt <- try(t.test(x, y, var.equal = FALSE), silent = TRUE)
      if (inherits(tt, "try-error")) NA_real_ else tt$p.value
    }, numeric(1))

    # BH correction across metabolites for this treatment-vs-control comparison
    padj <- p.adjust(pvals, method = "BH")

    # --- 4) Gate on the SAME effect shown in the heatmap: |LFC_heatmap| >= lfc_gate ---
    lfc_vec <- lfc[, smpl, drop = TRUE]  # log2((mean_trt+eps)/(mean_ctrl+eps))

    sig <- !is.na(padj) & (padj < alpha) & (abs(lfc_vec) >= lfc_gate)

    # Fill stars for this treatment column
    stars[, smpl] <- ifelse(sig, "*", "")
  }

  list(lfc = lfc, stars = stars)
}


#' Plot a faceted panel of metabolite log2 fold-changes vs control
#'
#' Creates a multi-panel (faceted) boxplot + jitter plot showing per-replicate
#' \eqn{\log_2} fold-changes (sample / CTRL) for a selected set of metabolites.
#' The function expects a wide data frame where replicate measurements are stored
#' in separate columns with names like \code{PREFIX_1}, \code{PREFIX_2}, etc.
#'
#' Workflow overview:
#' \itemize{
#'   \item Identify replicate columns using \code{replicate_regex}.
#'   \item Convert to long format and extract sample prefixes.
#'   \item Compute per-metabolite control means from \code{ctrl_prefix}.
#'   \item Normalize each replicate value to its metabolite's control mean and
#'     compute \eqn{\log_2} fold-change.
#'   \item Drop CTRL from the plotted samples (CTRL is the denominator), then plot
#'     fold-changes as faceted boxplots/jitter points.
#' }
#'
#' @param df Data frame containing replicate sample columns. Metabolites are
#'   identified by \code{rownames(df)}.
#' @param metabolites Character vector of metabolite names to plot. Only those
#'   present in \code{rownames(df)} are used.
#' @param ctrl_prefix Character. Control prefix (e.g., \code{"CTRL"}) used as the
#'   denominator for fold-change computations.
#' @param n_rows Integer. Number of facet rows.
#' @param n_cols Integer. Number of facet columns.
#' @param replicate_regex Character. Regex used to detect replicate columns.
#'   Default \code{"^[^_]+_\\\\d+$"} matches names like \code{"A_1"}.
#' @param tiny_pseudocount Numeric. Small value added to both numerator and
#'   denominator when computing relative abundance:
#'   \itemize{
#'     \item If \code{> 0}, uses \code{(Abundance + tiny_pseudocount) / (ctrl_mean + tiny_pseudocount)}
#'       to avoid division by zero and stabilize very small values.
#'     \item If \code{0}, uses \code{Abundance / ctrl_mean} and requires \code{ctrl_mean != 0}.
#'   }
#' @param y_limits Numeric length-2 vector. If not \code{NULL}, sets y-axis limits
#'   via \code{coord_cartesian(ylim = y_limits)}.
#' @param show_guides Logical. If \code{TRUE}, adds horizontal guide lines at
#'   \code{0} (solid) and \code{-2, 2} (dashed), useful as effect-size references.
#' @param palette Optional. Either:
#'   \itemize{
#'     \item unnamed vector of colors (must have at least as many colors as samples), or
#'     \item named vector mapping sample prefixes to colors.
#'   }
#'   Applied to both fill and color scales.
#' @param facet_label_width Integer. Wrap width for facet strip labels
#'   (metabolite names) using \code{stringr::str_wrap()}.
#' @param debug Logical. Reserved for debugging; currently unused.
#'
#' @details
#' **Metabolite matching**
#' Metabolites are pulled from \code{rownames(df)}. Requested metabolites not
#' present are ignored; if none are present, the function errors.
#'
#' **Control mean validity**
#' Metabolites are excluded from plotting if the control mean is not finite
#' (e.g., all NA) or if \code{tiny_pseudocount == 0} and \code{ctrl_mean == 0}
#' (to prevent division by zero). If all requested metabolites are excluded,
#' the function errors.
#'
#' **Plot semantics**
#' The plot shows distributions of replicate-level fold-changes per sample prefix.
#' CTRL replicates are used only to compute denominators and are not displayed.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_metabolites_lfc_panel(
#'   df = metab_wide,
#'   metabolites = c("Lactic acid", "Acetate", "Succinate"),
#'   ctrl_prefix = "CTRL",
#'   n_rows = 2, n_cols = 2,
#'   tiny_pseudocount = 1e-8,
#'   y_limits = c(-6, 6),
#'   palette = c(CTRL = "grey80", A = "#1b9e77", B = "#d95f02") # named palette example
#' )
#' p
#' }
#'
#' @export
plot_metabolites_lfc_panel <- function(df,
                                       metabolites,
                                       ctrl_prefix = "CTRL",
                                       n_rows = 2,
                                       n_cols = 3,
                                       replicate_regex = "^[^_]+_\\d+$",
                                       tiny_pseudocount = 0,
                                       y_limits = c(-10, 10),
                                       show_guides = TRUE,
                                       palette = NULL,
                                       facet_label_width = 28,
                                       debug = FALSE) {
  # --- Identify replicate/sample columns ---
  # Replicate columns are detected by regex (e.g., PREFIX_1, PREFIX_2, ...).
  sample_cols <- grep(replicate_regex, colnames(df), value = TRUE)
  if (!length(sample_cols)) stop("No replicate columns found.")

  # Derive sample prefixes from replicate column names
  unique_samples <- unique(sub("_[^_]+$", "", sample_cols))
  if (!ctrl_prefix %in% unique_samples) {
    stop(sprintf(
      "CTRL prefix '%s' not among samples: %s",
      ctrl_prefix, paste(unique_samples, collapse = ", ")
    ))
  }
  # Put control first, followed by other samples
  samp_levels <- c(ctrl_prefix, setdiff(unique_samples, ctrl_prefix))

  # --- Lift rownames & subset metabolites ---
  # Use rownames(df) as metabolite IDs and subset to requested metabolites.
  df2 <- as.data.frame(df[, sample_cols, drop = FALSE], stringsAsFactors = FALSE)
  df2$Metabolite <- rownames(df)

  metabolites <- as.character(metabolites)
  present <- intersect(metabolites, df2$Metabolite)
  if (!length(present)) stop("None of the requested metabolites found.")

  met_order <- present
  df2 <- dplyr::filter(df2, Metabolite %in% present)

  # --- Convert to long format ---
  # Long format has one row per Metabolite x Replicate column measurement.
  long <- df2 |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(sample_cols),
      names_to = "Replicate",
      values_to = "Abundance"
    ) |>
    dplyr::mutate(
      Sample     = sub("_[^_]+$", "", Replicate),
      Sample     = factor(Sample, levels = samp_levels),
      Metabolite = as.character(Metabolite),
      Abundance  = suppressWarnings(as.numeric(Abundance))
    )

  # --- Compute CTRL means per metabolite ---
  ctrl_means <- long |>
    dplyr::filter(Sample == ctrl_prefix) |>
    dplyr::group_by(Metabolite) |>
    dplyr::summarise(ctrl_mean = mean(Abundance, na.rm = TRUE), .groups = "drop")

  # Identify metabolites with invalid control means (non-finite, or zero without pseudocount)
  bad <- ctrl_means |>
    dplyr::filter(!is.finite(ctrl_mean) | (tiny_pseudocount == 0 & ctrl_mean == 0)) |>
    dplyr::pull(Metabolite)

  keep_mets <- setdiff(met_order, bad)
  if (!length(keep_mets)) stop("All requested metabolites had invalid CTRL means; nothing to plot.")

  # --- Normalize to CTRL and compute log2FC ---
  # Compute Relative = (Abundance + pc)/(ctrl_mean + pc) if pc>0, else Abundance/ctrl_mean.
  long2 <- long |>
    dplyr::filter(Metabolite %in% keep_mets) |>
    dplyr::left_join(ctrl_means, by = "Metabolite") |>
    dplyr::mutate(
      Relative = if (tiny_pseudocount > 0)
        (Abundance + tiny_pseudocount) / (ctrl_mean + tiny_pseudocount)
      else
        Abundance / ctrl_mean,
      log2FC = log2(Relative)
    ) |>
    dplyr::filter(is.finite(log2FC))

  # Remove CTRL from the plotted data (CTRL is the denominator, not a plotted group)
  long2 <- dplyr::filter(long2, Sample != ctrl_prefix)
  samp_levels <- setdiff(samp_levels, ctrl_prefix)
  long2 <- dplyr::mutate(long2, Metabolite = factor(Metabolite, levels = met_order))

  # --- Plot: boxplots + jitter, faceted by metabolite ---
  p <- ggplot2::ggplot(long2, ggplot2::aes(x = Sample, y = log2FC, fill = Sample)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = Sample),
      width = 0.2, size = 1.8, alpha = 0.85,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Metabolite,
      nrow = n_rows, ncol = n_cols,
      scales = "fixed", drop = FALSE,
      labeller = ggplot2::labeller(
        Metabolite = function(x) stringr::str_wrap(x, width = facet_label_width)
      )
    ) +
    ggplot2::ggtitle(paste("log2 Fold-Change vs", ctrl_prefix)) +
    ggplot2::labs(x = NULL, y = expression(log[2]("sample / CTRL"))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks.x = ggplot2::element_line(),
      legend.title = ggplot2::element_text(),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(title = "Sample"),
      color = "none"
    )

  # Optional y-axis limits (visual clipping only)
  if (!is.null(y_limits)) {
    p <- p +
      ggplot2::coord_cartesian(ylim = y_limits) +
      ggplot2::scale_y_continuous(breaks = pretty(y_limits, n = 7))
  }

  # Optional horizontal guide lines at 0 and +/-2
  if (show_guides) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "solid") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed")
  }

  # Optional palette (named or unnamed)
  if (!is.null(palette)) {
    if (is.null(names(palette))) {
      # Unnamed palette: assign in order of samp_levels
      if (length(palette) < length(samp_levels)) {
        stop(sprintf(
          "Palette has %d colors but needs at least %d.",
          length(palette), length(samp_levels)
        ))
      }
      pal_vec <- setNames(palette[seq_along(samp_levels)], samp_levels)
    } else {
      # Named palette: must cover all samples in samp_levels
      missing_cols <- setdiff(samp_levels, names(palette))
      if (length(missing_cols)) {
        stop(sprintf("Palette missing colors for: %s", paste(missing_cols, collapse = ", ")))
      }
      pal_vec <- palette[samp_levels]
    }

    p <- p +
      ggplot2::scale_fill_manual(values = pal_vec) +
      ggplot2::scale_color_manual(values = pal_vec)
  }

  return(p)
}



# ----- limma markers analysis -----
#' Align a metabolite feature table to sample metadata
#'
#' Internal helper to subset and reorder a metabolite matrix/data frame so that
#' its sample columns match the sample ordering in the metadata. This is a common
#' prerequisite for ordination, regression, differential testing, and plotting.
#'
#' The function:
#' \itemize{
#'   \item locates sample IDs in \code{metadata_df} (either in \code{sample_id_col},
#'     a \code{Sample} column, or rownames)
#'   \item finds overlapping sample IDs between \code{colnames(metab_df)} and metadata
#'   \item subsets \code{metab_df} to shared samples and converts it to a matrix
#'   \item reorders \code{metadata_df} rows to match the metabolite matrix column order
#' }
#'
#' @param metab_df Metabolite feature table with metabolites as rows and samples as
#'   columns. Column names must be sample IDs.
#' @param metadata_df Data frame containing sample metadata. Sample IDs are expected
#'   either in \code{sample_id_col}, in a column named \code{Sample}, or as rownames.
#' @param sample_id_col Optional. Column name in \code{metadata_df} containing sample IDs.
#'   If \code{NULL}, the function will use \code{Sample} if present; otherwise it will
#'   fall back to rownames and create a \code{Sample} column.
#'
#' @details
#' If \code{sample_id_col} is \code{NULL} and \code{metadata_df} has no \code{Sample}
#' column, rownames(\code{metadata_df}) are used as sample IDs and copied into a new
#' \code{Sample} column. This makes downstream code more consistent.
#'
#' Only samples present in both \code{metab_df} and \code{metadata_df} are retained.
#' The returned metadata is ordered to match the returned matrix columns.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{X}: numeric matrix of metabolites x aligned samples.
#'   \item \code{meta}: metadata data.frame reordered to match \code{colnames(X)}.
#' }
#'
#' @examples
#' \dontrun{
#' out <- .align_to_metadata(metab_df = mat, metadata_df = md, sample_id_col = "Sample")
#' X  <- out$X
#' md <- out$meta
#' stopifnot(identical(colnames(X), rownames(md)))
#' }
#'
#' @keywords internal
.align_to_metadata <- function(metab_df, metadata_df, sample_id_col = NULL) {
  # Require sample IDs in the feature table column names
  stopifnot(!is.null(colnames(metab_df)))
  md <- metadata_df

  # Determine where sample IDs live in metadata
  if (is.null(sample_id_col)) {
    if ("Sample" %in% colnames(md)) {
      sample_id_col <- "Sample"
    } else if (!is.null(rownames(md))) {
      md$Sample <- rownames(md)
      sample_id_col <- "Sample"
    } else {
      stop("metadata_df must have a 'Sample' column or rownames = sample IDs.")
    }
  }

  # Identify overlapping samples
  common <- intersect(colnames(metab_df), md[[sample_id_col]])
  if (!length(common)) stop("No overlapping sample IDs between metabolome and metadata.")

  # Subset metabolite table to shared samples and convert to matrix
  X <- as.matrix(metab_df[, common, drop = FALSE])

  # Reorder metadata to match the metabolite matrix columns
  md <- md[match(common, md[[sample_id_col]]), , drop = FALSE]
  rownames(md) <- md[[sample_id_col]]

  list(X = X, meta = md)
}



#' Identify metabolite markers by cluster using limma (one-vs-rest and optional pairwise)
#'
#' Fits a linear model per metabolite using \code{limma} to test for differential
#' abundance across clusters, optionally adjusting for covariates and optionally
#' accounting for within-block correlation (e.g., repeated measures / SynCom blocks)
#' using \code{limma::duplicateCorrelation()}.
#'
#' The function returns:
#' \itemize{
#'   \item one-vs-rest contrasts for each cluster (each cluster vs the average of all others),
#'   \item optional pairwise contrasts among all cluster levels,
#'   \item the design matrix used and the estimated duplicate correlation (if any).
#' }
#'
#' @param metab_df Metabolite feature table with metabolites as rows and samples as
#'   columns. Column names must be sample IDs.
#' @param metadata_df Sample metadata data frame. Sample IDs are aligned to
#'   \code{metab_df} using \code{.align_to_metadata()}.
#' @param sample_id_col Optional. Column in \code{metadata_df} containing sample IDs.
#'   If \code{NULL}, \code{.align_to_metadata()} will try \code{"Sample"} or rownames.
#' @param cluster_var Character. Metadata column defining clusters (e.g.,
#'   \code{"ATTRIBUTE_Cluster"}). Treated as a factor.
#' @param covariates Optional character vector of metadata columns to include as
#'   covariates (e.g., \code{c("ATTRIBUTE_Time")}). Character columns are coerced
#'   to factors; numeric columns are kept numeric. Dummy coding is generated via
#'   \code{model.matrix(~ 0 + ., data = cov_df)}.
#' @param block_var Optional character. Metadata column defining a blocking factor
#'   for duplicate correlation (e.g., \code{"ATTRIBUTE_SynCom"}). If provided,
#'   \code{duplicateCorrelation} is used and the model is fit with \code{block}.
#' @param log_transform Logical. If \code{TRUE}, apply \code{log(X + log_offset)} to
#'   the aligned data prior to modeling (default \code{TRUE}). If data contain
#'   negative values, they are shifted to be non-negative before logging.
#' @param log_offset Numeric. Offset added before log transform (default \code{1}).
#' @param do_pairwise Logical. If \code{TRUE}, also compute all pairwise contrasts
#'   among cluster levels (default \code{TRUE}).
#' @param adjust_method Multiple testing adjustment passed to \code{limma::topTable()}
#'   (default \code{"BH"}).
#'
#' @details
#' **Design**
#' The cluster effect is encoded as a no-intercept design:
#' \code{~ 0 + cluster}, producing one coefficient per cluster level. If covariates
#' are provided, their dummy-coded columns are appended to the design matrix.
#'
#' **One-vs-rest contrasts**
#' For \eqn{K} clusters, each contrast is:
#' \deqn{c_i = 1 \text{ for cluster } i,\;\; -1/(K-1) \text{ for all other clusters}}
#' so that each cluster is compared to the average of the rest.
#'
#' **Pairwise contrasts**
#' If enabled, all combinations of cluster levels are tested as \code{B - A}.
#'
#' **Blocking / correlation**
#' If \code{block_var} is provided, \code{duplicateCorrelation()} estimates a
#' consensus correlation to account for non-independence within blocks, and that
#' correlation is used in \code{lmFit()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{markers_one_vs_rest}: data frame of limma results for one-vs-rest
#'     contrasts. Includes \code{Metabolite}, \code{Contrast}, and \code{TargetCluster}.
#'   \item \code{markers_pairwise}: data frame of limma results for pairwise contrasts,
#'     or \code{NULL} if disabled.
#'   \item \code{duplicate_correlation}: estimated consensus correlation (or \code{NA}
#'     if no blocking was used).
#'   \item \code{design}: the final design matrix used for modeling.
#'   \item \code{cluster_levels}: character vector of cluster level names.
#'   \item \code{used_args}: list of key arguments used for reproducibility.
#' }
#'
#' @examples
#' \dontrun{
#' res <- limma_markers_by_cluster_general(
#'   metab_df = metab_mat,
#'   metadata_df = md,
#'   cluster_var = "ATTRIBUTE_Cluster",
#'   covariates = c("ATTRIBUTE_Time"),
#'   block_var = "ATTRIBUTE_SynCom",
#'   log_transform = TRUE,
#'   log_offset = 1,
#'   do_pairwise = TRUE
#' )
#' head(res$markers_one_vs_rest)
#' head(res$markers_pairwise)
#' }
#'
#' @export
limma_markers_by_cluster_general <- function(
    metab_df, metadata_df,
    sample_id_col = NULL,          # e.g., "Sample" (or set rownames(metadata) accordingly and pass "Sample")
    cluster_var,            # e.g., "ATTRIBUTE_Cluster"
    covariates = NULL,      # e.g., c("ATTRIBUTE_Time")
    block_var = NULL,       # e.g., "ATTRIBUTE_SynCom"  (optional)
    log_transform = TRUE, log_offset = 1,
    do_pairwise = TRUE,
    adjust_method = "BH"
) {
  # ---- Dependency check ----
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Install it via install.packages('limma').")
  }

  # ---- Align metabolite table to metadata ----
  aligned <- .align_to_metadata(metab_df, metadata_df, sample_id_col)
  X  <- aligned$X
  md <- aligned$meta

  # ---- Build model pieces ----
  if (!cluster_var %in% colnames(md)) stop("cluster_var not found in metadata.")
  cluster <- factor(md[[cluster_var]])

  # Cluster design: one coefficient per cluster (no intercept)
  Z_cluster <- model.matrix(~ 0 + cluster)
  colnames(Z_cluster) <- paste0("cluster", levels(cluster))  # e.g., cluster1, cluster2...

  # Covariate design (optional): character -> factor; numeric kept numeric
  cov_mat <- NULL
  if (!is.null(covariates) && length(covariates) > 0) {
    miss <- setdiff(covariates, colnames(md))
    if (length(miss)) stop("Missing covariate(s) in metadata: ", paste(miss, collapse = ", "))

    cov_df <- md[, covariates, drop = TRUE]
    cov_df <- as.data.frame(cov_df, stringsAsFactors = FALSE)

    for (nm in names(cov_df)) {
      if (is.character(cov_df[[nm]])) cov_df[[nm]] <- factor(cov_df[[nm]])
    }

    # Dummy-code factors and keep numeric predictors
    cov_mat <- model.matrix(~ 0 + ., data = cov_df)
    colnames(cov_mat) <- paste0("cov_", colnames(cov_mat))
  }

  # Combine cluster and covariates into final design
  design <- if (is.null(cov_mat)) Z_cluster else cbind(Z_cluster, cov_mat)

  # ---- Optional log transform ----
  # Shifts negative values to be non-negative before log transform.
  if (log_transform) {
    minX <- suppressWarnings(min(X, na.rm = TRUE))
    if (is.finite(minX) && minX < 0) X <- X - minX
    X <- log(X + log_offset)
  }

  # ---- Fit limma model (optionally with blocking) ----
  use_block <- !is.null(block_var)
  if (use_block) {
    if (!block_var %in% colnames(md)) stop("block_var not found in metadata.")
    block <- md[[block_var]]

    # Estimate within-block correlation and fit model with block/correlation
    corfit <- limma::duplicateCorrelation(X, design, block = block)
    fit <- limma::lmFit(X, design, block = block, correlation = corfit$consensus)
  } else {
    corfit <- list(consensus = NA_real_)
    fit <- limma::lmFit(X, design)
  }

  # ---- One-vs-rest contrasts over cluster columns ----
  cl_cols <- colnames(Z_cluster)   # e.g., cluster1, cluster2, ...
  K <- length(cl_cols)
  if (K < 2) stop("Need at least 2 cluster levels.")

  # For each cluster i: 1 for i, -1/(K-1) for others
  L_ovr <- sapply(seq_len(K), function(i) {
    v <- rep(-1/(K - 1), K)
    v[i] <- 1
    names(v) <- cl_cols
    v
  })
  colnames(L_ovr) <- paste0(sub("^cluster", "", cl_cols), "_vs_rest")

  # Embed contrasts into the full design space (zeros for covariate coefficients)
  L_full <- matrix(
    0,
    nrow = ncol(design),
    ncol = ncol(L_ovr),
    dimnames = list(colnames(design), colnames(L_ovr))
  )
  L_full[cl_cols, colnames(L_ovr)] <- L_ovr

  # Fit contrasts and apply empirical Bayes moderation
  fit_ovr <- limma::contrasts.fit(fit, L_full)
  fit_ovr <- limma::eBayes(fit_ovr)

  # Collect results for all one-vs-rest contrasts
  markers_ovr <- lapply(seq_len(ncol(L_full)), function(i) {
    limma::topTable(fit_ovr, coef = i, number = Inf, adjust.method = adjust_method) %>%
      rownames_to_column("Metabolite") %>%
      mutate(
        Contrast = colnames(L_full)[i],
        TargetCluster = sub("_vs_rest$", "", colnames(L_full)[i])
      )
  }) %>% bind_rows()

  # ---- Optional: pairwise contrasts among clusters ----
  markers_pairwise <- NULL
  if (isTRUE(do_pairwise) && K >= 2) {
    combs <- t(combn(cl_cols, 2))

    # Build B - A contrasts for every pair (A,B)
    L_pw <- sapply(seq_len(nrow(combs)), function(i) {
      v <- setNames(rep(0, K), cl_cols)
      A <- combs[i, 1]
      B <- combs[i, 2]
      v[B] <-  1
      v[A] <- -1
      v
    })
    colnames(L_pw) <- apply(combs, 1, function(x) {
      paste0(sub("^cluster", "", x[2]), "_vs_", sub("^cluster", "", x[1]))
    })

    # Embed pairwise contrasts into full design space
    Lpw_full <- matrix(
      0,
      nrow = ncol(design),
      ncol = ncol(L_pw),
      dimnames = list(colnames(design), colnames(L_pw))
    )
    Lpw_full[cl_cols, colnames(L_pw)] <- L_pw

    # Fit pairwise contrasts and apply empirical Bayes moderation
    fit_pw <- limma::contrasts.fit(fit, Lpw_full)
    fit_pw <- eBayes(fit_pw)

    # Collect results for all pairwise contrasts
    markers_pairwise <- lapply(seq_len(ncol(Lpw_full)), function(i) {
      topTable(fit_pw, coef = i, number = Inf, adjust.method = adjust_method) %>%
        rownames_to_column("Metabolite") %>%
        mutate(Contrast = colnames(Lpw_full)[i])
    }) %>% bind_rows()
  }

  list(
    markers_one_vs_rest     = markers_ovr,
    markers_pairwise        = markers_pairwise,
    duplicate_correlation   = corfit$consensus,
    design                  = design,
    cluster_levels          = sub("^cluster", "", cl_cols),
    used_args               = list(
      cluster_var = cluster_var,
      covariates = covariates,
      block_var = block_var,
      log_transform = log_transform,
      log_offset = log_offset,
      adjust_method = adjust_method
    )
  )
}


#' Summarize limma cluster markers and draw a heatmap annotated with SIRIUS classes
#'
#' Selects top-N marker metabolites per cluster from a limma one-vs-rest result,
#' aligns the metabolite matrix to metadata, optionally log-transforms and row-scales
#' the selected metabolite profiles, and renders a \code{ComplexHeatmap} heatmap
#' with:
#' \itemize{
#'   \item column annotations for sample cluster membership, and
#'   \item row annotations derived from SIRIUS/ClassyFire chemical class columns.
#' }
#'
#' Optionally saves the heatmap to a file (PDF/SVG/PNG/TIFF/JPG) with configurable
#' dimensions, DPI, legend placement, and legend merging.
#'
#' @param metab_df Metabolite feature table with metabolites as rows and samples as
#'   columns. Column names must be sample IDs.
#' @param metadata_df Sample metadata data frame containing \code{cluster_var} and
#'   sample identifiers (rownames or \code{sample_id_col}).
#' @param sample_id_col Optional. Column name in \code{metadata_df} holding sample IDs.
#'   If \code{NULL}, non-empty rownames are used.
#' @param cluster_var Character. Metadata column defining sample clusters (e.g.,
#'   \code{"ATTRIBUTE_Cluster"}).
#' @param sirius_df Data frame of SIRIUS annotations, including an identifier column
#'   (see \code{id_col}) and one or more class columns (see \code{class_cols}).
#' @param id_col Character. Column name in \code{sirius_df} containing the identifier
#'   used to match metabolites to SIRIUS rows (default \code{"row.ID"}).
#' @param class_cols Character vector of SIRIUS/ClassyFire columns to use as row
#'   annotations. Only columns present in \code{sirius_df} are used.
#' @param id_pattern Regular expression used to extract an ID from metabolite rownames
#'   of \code{metab_df}. The default \code{"^X(\\\\d+).*"} extracts the first integer
#'   after an \code{"X"} (e.g., \code{"X1234_something"} \eqn{\rightarrow} \code{"1234"}).
#' @param limma_res Result list from \code{limma_markers_by_cluster_general()}, which must
#'   contain \code{$markers_one_vs_rest}.
#' @param top_n Integer. Number of top markers to select per cluster after filtering.
#' @param p_adj_thresh Numeric. Adjusted p-value cutoff (\code{adj.P.Val}) used to filter
#'   markers before selecting \code{top_n}.
#' @param min_logFC Numeric. Minimum \code{logFC} required to retain a marker (default \code{0}).
#' @param log_transform Logical. If \code{TRUE}, apply \code{log(X + log_offset)} to the
#'   selected matrix before plotting. Negative values (if any) are shifted to be non-negative.
#' @param log_offset Numeric. Offset used in the log transform (default \code{1}).
#' @param scale_rows Logical. If \code{TRUE}, row-scale selected metabolites for display
#'   (z-scores across samples). Non-finite z-scores are set to 0.
#' @param heatmap_colors Deprecated/unused in the current implementation (color mapping is
#'   defined via \code{circlize::colorRamp2()}).
#' @param cluster_colors Optional named vector mapping cluster levels to colors for the
#'   column annotation. If \code{NULL}, a \code{RColorBrewer::Set2} palette is used.
#' @param class_na_label Character. Placeholder label used when SIRIUS class annotations are
#'   missing/blank (default \code{"Unclassified"}).
#' @param class_na_color Color used for \code{class_na_label} (default \code{"#BDBDBD"}).
#' @param out_file Optional. If provided, saves the heatmap to this file path. Supported
#'   extensions: \code{pdf}, \code{svg}, \code{png}, \code{tif/tiff}, \code{jpg/jpeg}.
#' @param out_width Numeric. Output width in inches.
#' @param out_height Numeric. Output height in inches.
#' @param out_dpi Integer. DPI for raster outputs (PNG/JPG/TIFF).
#' @param c_legend_ncol Integer. Number of legend columns for column annotations.
#' @param r_legend_ncol Integer. Number of legend columns for row annotations.
#' @param legend_side Character. Legend placement for \code{ComplexHeatmap::draw()}:
#'   \code{"bottom"}, \code{"top"}, \code{"left"}, or \code{"right"}.
#' @param merge_legends Logical. If \code{TRUE}, merge heatmap and annotation legends into a
#'   single legend block (passed to \code{draw(merge_legend = ...)}).
#'
#' @details
#' **Step A: Sample alignment**
#' Samples are intersected between \code{colnames(metab_df)} and metadata IDs. The
#' metabolite matrix and metadata are subset and ordered consistently, and a column
#' annotation data frame is created with \code{Cluster = factor(metadata_df[[cluster_var]])}.
#'
#' **Step B: SIRIUS row annotations**
#' Metabolite rownames are mapped to SIRIUS IDs using \code{id_pattern}. Only distinct
#' SIRIUS rows per ID are kept. Selected class columns are joined onto metabolites.
#' Missing/blank class values are replaced by \code{class_na_label}.
#'
#' **Step C: Marker selection**
#' Uses \code{limma_res$markers_one_vs_rest} and filters by \code{adj.P.Val <= p_adj_thresh}
#' and \code{logFC >= min_logFC}. For each target cluster, the top \code{top_n} metabolites
#' are selected by increasing adjusted p-value and decreasing logFC.
#'
#' **Step D: Heatmap**
#' Columns are ordered by cluster. If \code{scale_rows = TRUE}, z-scores are computed per
#' metabolite across samples and used for the heatmap values. A symmetric blue-white-red
#' color mapping is built around 0. Column and row annotations are added with independent
#' legend column controls.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{heatmap}: \code{ComplexHeatmap::Heatmap} object.
#'   \item \code{X_display}: matrix used for plotting (optionally row-scaled).
#'   \item \code{X_values}: selected metabolite matrix values before scaling.
#'   \item \code{ann_col}: column annotation data.frame (clusters).
#'   \item \code{ann_row}: row annotation data.frame (SIRIUS classes).
#'   \item \code{annotation_colors}: list of palettes used for annotations (for reference).
#'   \item \code{top_table}: data.frame of selected top markers per cluster.
#'   \item \code{selected_metabolites}: character vector of metabolite IDs used in the heatmap.
#' }
#'
#' @examples
#' \dontrun{
#' res <- summarize_markers_and_heatmap_with_classes(
#'   metab_df = metab_mat,
#'   metadata_df = md,
#'   sample_id_col = "Sample",
#'   cluster_var = "ATTRIBUTE_Cluster",
#'   sirius_df = sirius_tbl,
#'   limma_res = limma_out,
#'   top_n = 10,
#'   p_adj_thresh = 0.05,
#'   scale_rows = TRUE,
#'   out_file = "markers_heatmap.pdf"
#' )
#' res$heatmap
#' res$top_table
#' }
#'
#' @export
summarize_markers_and_heatmap_with_classes <- function(
    # --- core inputs ---
  metab_df, metadata_df,
  sample_id_col = NULL,             # e.g. "Sample"; if NULL, uses metadata rownames
  cluster_var,                       # e.g. "ATTRIBUTE_Cluster"
  # --- SIRIUS inputs ---
  sirius_df,
  id_col = "row.ID",
  class_cols = c("SIRIUS_ClassyFire.class",
                 "SIRIUS_ClassyFire.most.specific.class",
                 "SIRIUS_ClassyFire.subclass",
                 "SIRIUS_ClassyFire.level.5"),
  id_pattern = "^X(\\d+).*",
  # --- limma selection ---
  limma_res,                         # output of limma_markers_by_cluster_general()
  top_n = 15,
  p_adj_thresh = 0.05,
  min_logFC = 0,
  log_transform = TRUE, log_offset = 1,
  scale_rows = TRUE,
  heatmap_colors = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),
  # --- plotting options ---
  cluster_colors = NULL,             # named vector for sample clusters (optional)
  class_na_label = "Unclassified",
  class_na_color = "#BDBDBD",
  out_file   = NULL,  # e.g. "markers_heatmap.pdf" or "markers_heatmap.png"
  out_width  = 9,     # inches
  out_height = 12,    # inches
  out_dpi    = 300,   # used for raster formats (png/jpg/tiff)
  c_legend_ncol = 2,     # columns for both heatmap & annotation legends,
  r_legend_ncol = 2,
  legend_side = "bottom",   # "bottom", "top", "left", or "right"
  merge_legends = FALSE     # TRUE = put all legends into one combined block
) {
  ## ---------- Step A: align + column annotation ----------
  # Validate required structures/columns and align samples between metab_df and metadata_df.
  if (is.null(colnames(metab_df))) stop("metab_df must have sample column names.")
  if (!(cluster_var %in% colnames(metadata_df))) {
    stop("Cluster column '", cluster_var, "' not found in metadata_df.")
  }

  # Determine where sample IDs live in metadata (explicit column preferred, else rownames)
  if (!is.null(sample_id_col) && sample_id_col %in% colnames(metadata_df)) {
    md_ids <- as.character(metadata_df[[sample_id_col]])
    md_id_src <- paste0("column '", sample_id_col, "'")
  } else if (!is.null(rownames(metadata_df)) && all(nchar(rownames(metadata_df)) > 0)) {
    md_ids <- rownames(metadata_df)
    md_id_src <- "metadata rownames"
    sample_id_col <- "<rownames>"
  } else {
    stop("Could not find sample IDs in metadata_df (need sample_id_col or non-empty rownames).")
  }

  # Restrict to overlapping samples
  common <- intersect(colnames(metab_df), md_ids)
  if (length(common) == 0) {
    stop("No overlapping sample IDs between metab_df and metadata_df (", md_id_src, ").")
  }

  # Subset metabolite table to common samples and align metadata row order
  X <- as.matrix(metab_df[, common, drop = FALSE])
  if (sample_id_col == "<rownames>") {
    md <- metadata_df[match(common, rownames(metadata_df)), , drop = FALSE]
    rownames(md) <- rownames(md)
  } else {
    md <- metadata_df[match(common, metadata_df[[sample_id_col]]), , drop = FALSE]
    rownames(md) <- md[[sample_id_col]]
  }

  # Column annotation: cluster membership per sample
  ann_col_all <- data.frame(Cluster = factor(md[[cluster_var]]))
  rownames(ann_col_all) <- rownames(md)
  message(
    "Prepared alignment: ", ncol(X), " samples; clusters = {",
    paste(levels(ann_col_all$Cluster), collapse = ", "),
    "}"
  )

  ## ---------- Step B: SIRIUS row annotations ----------
  # Validate requested SIRIUS ID and class columns; join class info to metabolites.
  if (!id_col %in% colnames(sirius_df)) {
    stop(
      "SIRIUS id_col '", id_col, "' not found in sirius_df. Available: ",
      paste(colnames(sirius_df), collapse = ", ")
    )
  }
  class_cols_present <- intersect(class_cols, colnames(sirius_df))
  if (length(class_cols_present) == 0) {
    stop(
      "None of the requested class_cols found in sirius_df. Available: ",
      paste(colnames(sirius_df), collapse = ", ")
    )
  }

  metas <- rownames(X)
  if (is.null(metas)) stop("metab_df must have metabolite rownames.")

  # Extract IDs from metabolite rownames for joining to sirius_df
  ids_extracted <- sub(id_pattern, "\\1", metas)

  # Keep one row per SIRIUS ID and only the requested class columns
  sirius_min <- sirius_df |>
    dplyr::mutate(.id = as.character(.data[[id_col]])) |>
    dplyr::distinct(.id, .keep_all = TRUE) |>
    dplyr::select(.id, dplyr::all_of(class_cols_present))

  # Join metabolite->ID mapping to SIRIUS classes
  id_map <- data.frame(Metabolite = metas, .id = ids_extracted, stringsAsFactors = FALSE)
  ann_row_full <- id_map |>
    dplyr::left_join(sirius_min, by = ".id") |>
    as.data.frame()
  rownames(ann_row_full) <- ann_row_full$Metabolite
  ann_row_full$Metabolite <- NULL
  if (".id" %in% colnames(ann_row_full)) ann_row_full$.id <- NULL

  # Drop annotation columns that are entirely NA; coerce remaining to factors
  if (ncol(ann_row_full) > 0) {
    all_na <- vapply(ann_row_full, function(x) all(is.na(x)), logical(1))
    if (any(all_na)) ann_row_full <- ann_row_full[, !all_na, drop = FALSE]
    for (nm in colnames(ann_row_full)) ann_row_full[[nm]] <- droplevels(factor(ann_row_full[[nm]]))
  }

  match_count <- sum(ids_extracted %in% unique(sirius_min$.id))
  message("SIRIUS matched IDs: ", match_count, " / ", length(ids_extracted))

  ## ---------- Step C: select top-N markers per cluster (from limma) ----------
  if (!("markers_one_vs_rest" %in% names(limma_res))) {
    stop("limma_res must contain $markers_one_vs_rest.")
  }

  ovr <- limma_res$markers_one_vs_rest
  req_cols <- c("Metabolite","adj.P.Val","logFC","TargetCluster")
  if (!all(req_cols %in% colnames(ovr))) {
    stop("limma_res$markers_one_vs_rest must have columns: ", paste(req_cols, collapse = ", "))
  }

  # Filter markers by thresholds, then take top_n per target cluster
  top_by_cluster <- ovr |>
    dplyr::filter(adj.P.Val <= p_adj_thresh, logFC >= min_logFC) |>
    dplyr::group_by(TargetCluster) |>
    dplyr::arrange(adj.P.Val, dplyr::desc(logFC), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  top_metabs <- unique(top_by_cluster$Metabolite)
  if (!length(top_metabs)) stop("No metabolites passed thresholds; relax p_adj_thresh or min_logFC.")

  # Subset feature table to selected metabolites
  keep <- intersect(rownames(X), top_metabs)
  if (!length(keep)) stop("Selected metabolites not found in rownames(metab_df). Check naming.")
  Xsub <- X[keep, , drop = FALSE]

  # Optional transform for display (match limma transform settings)
  if (isTRUE(log_transform)) {
    minX <- suppressWarnings(min(Xsub, na.rm = TRUE))
    if (is.finite(minX) && minX < 0) Xsub <- Xsub - minX
    Xsub <- log(Xsub + log_offset)
  }

  # Order columns by cluster level for visualization
  md[[cluster_var]] <- factor(md[[cluster_var]])
  col_order <- order(md[[cluster_var]], decreasing = FALSE)
  Xsub <- Xsub[, col_order, drop = FALSE]
  ann_col <- ann_col_all[colnames(Xsub), , drop = FALSE]

  # Subset row annotations to selected metabolites
  ann_row <- ann_row_full[rownames(Xsub), , drop = FALSE]

  # Replace missing/blank class values with a placeholder label and drop unused levels
  if (ncol(ann_row) > 0) {
    for (nm in colnames(ann_row)) {
      v <- as.character(ann_row[[nm]])
      v[is.na(v) | trimws(v) == ""] <- class_na_label
      ann_row[[nm]] <- droplevels(factor(v))
    }
  }

  # Optional row scaling for display (z-scores per metabolite)
  X_display <- Xsub
  if (isTRUE(scale_rows) && nrow(X_display) > 1) {
    X_display <- t(scale(t(X_display)))
    X_display[!is.finite(X_display)] <- 0
  }

  # ---------- Step D: build annotation colors & plot ----------
  # Column (Cluster) palette
  if (is.null(cluster_colors)) {
    k <- nlevels(ann_col$Cluster)
    cluster_colors <- setNames(
      brewer.pal(max(3, min(8, k)), "Set2")[seq_len(k)],
      levels(ann_col$Cluster)
    )
  }

  # Store palettes in a single list for reference / downstream reuse
  ann_colors <- list(Cluster = cluster_colors)

  # Row annotation palettes: discrete, include placeholder label color
  if (!is.null(ann_row) && ncol(ann_row) > 0) {
    make_discrete_pal <- function(vals) {
      lev <- levels(vals)

      has_na_lab <- class_na_label %in% lev
      lev_core   <- setdiff(lev, class_na_label)
      n <- length(lev_core)

      base <- if (n <= 12) RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)] else scales::hue_pal()(n)
      pal <- setNames(base, lev_core)
      if (has_na_lab) pal[[class_na_label]] <- class_na_color
      pal
    }
    row_palettes <- lapply(ann_row, make_discrete_pal)
    ann_colors   <- c(ann_colors, row_palettes)
  }

  # ---- ComplexHeatmap dependencies ----
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install packages 'ComplexHeatmap' and 'circlize' for multi-column legends.")
  }

  # Heatmap color function (symmetric around 0 for z-scored data)
  zlim <- max(abs(range(X_display, finite = TRUE)))
  col_fun <- circlize::colorRamp2(c(-zlim, 0, zlim), c("#2166AC","#F7F7F7","#B2182B"))

  # Column annotation colors
  col_ann_cols <- list(Cluster = cluster_colors)

  # Row annotation colors: use the palettes built above so placeholder label gets its color
  row_ann_cols <- if (!is.null(ann_row) && ncol(ann_row) > 0) ann_colors[names(ann_row)] else list()

  # Top (column) annotation with legend column control
  ha_top <- ComplexHeatmap::HeatmapAnnotation(
    df  = ann_col,
    col = col_ann_cols,
    annotation_legend_param = list(ncol = c_legend_ncol)
  )

  # Left (row) annotation with legend column control (or NULL if no row annotations)
  ha_left <- if (ncol(ann_row) > 0) {
    ComplexHeatmap::rowAnnotation(
      df  = ann_row,
      col = row_ann_cols,
      annotation_legend_param = list(ncol = r_legend_ncol)
    )
  } else {
    NULL
  }

  # Main heatmap object
  ht <- ComplexHeatmap::Heatmap(
    X_display,
    name = "z",
    col  = col_fun,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    top_annotation  = ha_top,
    left_annotation = ha_left,
    heatmap_legend_param = list(ncol = 2),
    column_title = sprintf("Top markers per cluster (top %d, FDR <= %.02f)", top_n, p_adj_thresh)
  )

  # ---- Helper for saving heatmaps to common formats ----
  save_ht <- function(ht, file, width, height, dpi, legend_side, merge_legends) {
    ext <- tolower(tools::file_ext(file))

    if (ext %in% c("pdf")) {
      grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    } else if (ext %in% c("svg")) {
      grDevices::svg(file, width = width, height = height)
    } else if (ext %in% c("png")) {
      grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi, units = "px")
    } else if (ext %in% c("tif","tiff")) {
      grDevices::tiff(file, width = width * dpi, height = height * dpi, res = dpi, units = "px", compression = "lzw")
    } else if (ext %in% c("jpg","jpeg")) {
      grDevices::jpeg(file, width = width * dpi, height = height * dpi, res = dpi, units = "px", quality = 95)
    } else {
      stop("Unsupported file extension: ", ext)
    }
    on.exit(grDevices::dev.off(), add = TRUE)

    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = legend_side,
      annotation_legend_side = legend_side,
      merge_legend = merge_legends
    )
  }

  # ---- Draw or save ----
  if (!is.null(out_file)) {
    save_ht(ht, out_file, out_width, out_height, out_dpi, legend_side, merge_legends)
  } else {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = legend_side,
      annotation_legend_side = legend_side,
      merge_legend = merge_legends
    )
  }

  # Return key objects for reuse/inspection
  list(
    heatmap              = ht,
    X_display            = X_display,
    X_values             = Xsub,
    ann_col              = ann_col,
    ann_row              = ann_row,
    annotation_colors    = ann_colors,
    top_table            = top_by_cluster,
    selected_metabolites = rownames(Xsub)
  )
}


#' Load a BIOM file as a taxonomically agglomerated feature table
#'
#' Imports a BIOM file (containing an OTU/ASV table and taxonomy), agglomerates
#' features to a requested Greengenes-style taxonomic rank (e.g., Genus/Species),
#' and returns a clean numeric table suitable for downstream ecological analyses
#' (e.g., \pkg{vegan}) and visualization (barplots/heatmaps).
#'
#' The function:
#' \itemize{
#'   \item validates/derives the expected taxonomy column names for \code{tax_rank}
#'     (via \code{get_colNames_per_rank()})
#'   \item imports the BIOM using \code{phyloseq::import_biom()}
#'   \item parses taxonomy using either Greengenes defaults or a custom strain parser
#'   \item agglomerates/aggregates the feature table at the requested rank
#'     (via \code{extract_table()})
#'   \item cleans the resulting table (unique rownames, optional ordering)
#'     (via \code{clean_table()})
#' }
#'
#' @param biom_path Character. Path to a \code{.biom} file containing a feature table
#'   and taxonomy annotations (e.g., exported from QIIME2/HMP pipelines).
#' @param tax_rank Character. Taxonomic rank to aggregate to (default \code{"Species"}).
#'   Must be a valid Greengenes rank string expected by your helper functions
#'   (e.g., Kingdom/Phylum/Class/Order/Family/Genus/Species), otherwise downstream
#'   helpers may error.
#' @param strain_taxonomy Logical. If \code{TRUE}, taxonomy strings are parsed with
#'   \code{parse_taxonomy_strain} to preserve strain-level taxonomy conventions.
#'   If \code{FALSE}, uses \code{phyloseq::parse_taxonomy_greengenes}.
#' @param order_table Logical. If \code{TRUE}, order rows from larger to smaller
#'   mean abundance (typically via \code{rowMeans}) during cleaning.
#'
#' @details
#' **Output format**
#' The returned object is a data frame/matrix where:
#' \itemize{
#'   \item rows correspond to aggregated taxa at \code{tax_rank}
#'   \item columns correspond to samples
#'   \item rownames are the taxon labels at \code{tax_rank} (dereplicated so they are unique)
#' }
#'
#' The exact agglomeration and cleaning logic is delegated to project helper
#' functions \code{extract_table()} and \code{clean_table()}.
#'
#' @return A cleaned abundance table (data frame) with taxa as rows and samples as
#'   columns, ready for use in downstream analyses.
#'
#' @examples
#' \dontrun{
#' # Species-level table using Greengenes taxonomy parsing
#' tab_species <- load_biom_as_table(
#'   biom_path = "hmp_table.biom",
#'   tax_rank = "Species",
#'   strain_taxonomy = FALSE,
#'   order_table = TRUE
#' )
#'
#' # Genus-level table with custom strain taxonomy parsing
#' tab_genus <- load_biom_as_table(
#'   biom_path = "hmp_table.biom",
#'   tax_rank = "Genus",
#'   strain_taxonomy = TRUE
#' )
#' }
#'
#' @export
load_biom_as_table <- function(biom_path,
                               tax_rank = "Species",
                               strain_taxonomy = FALSE,
                               order_table = FALSE) {

  # Determine the taxonomy column(s) to use/construct for the requested rank
  unite_colNames <- get_colNames_per_rank(tax_rank)

  # Import BIOM with the appropriate taxonomy parsing function
  if (isTRUE(strain_taxonomy)) {
    biom_object <- phyloseq::import_biom(biom_path, parseFunction = parse_taxonomy_strain)
  } else {
    biom_object <- phyloseq::import_biom(
      biom_path,
      parseFunction = phyloseq::parse_taxonomy_greengenes
    )
  }

  # Agglomerate/reshape the feature table at the requested taxonomic rank
  extracted_feature_table <- extract_table(biom_object, tax_rank, unite_colNames)

  # Clean table (unique taxa names; optional ordering by mean abundance)
  return(clean_table(extracted_feature_table, order_table = order_table))
}


#' Get taxonomy column names to concatenate for a given Greengenes rank
#'
#' Returns the set of taxonomy table column names that should be joined/concatenated
#' to form a taxonomic label at the requested rank. This is used when agglomerating
#' OTU/ASV tables from BIOM/phyloseq objects to a higher taxonomic level.
#'
#' Supported ranks (as implemented here):
#' \itemize{
#'   \item \code{"Strain"}  \eqn{\rightarrow} \code{c("Genus","Species","Strain")}
#'   \item \code{"Species"} \eqn{\rightarrow} \code{c("Genus","Species")}
#'   \item \code{"Genus"}   \eqn{\rightarrow} \code{c("Genus")}
#'   \item \code{"Family"}  \eqn{\rightarrow} \code{c("Family")}
#'   \item \code{"Order"}   \eqn{\rightarrow} \code{c("Order")}
#' }
#'
#' @param tax_rank Character. Taxonomic rank in Greengenes-style naming.
#'
#' @return Character vector of taxonomy column names to use for building labels at
#'   the requested rank.
#'
#' @examples
#' \dontrun{
#' get_colNames_per_rank("Species")  # c("Genus","Species")
#' get_colNames_per_rank("Genus")    # c("Genus")
#' }
#'
#' @export
get_colNames_per_rank <- function(tax_rank) {
  # Map rank -> taxonomy columns required to construct a label at that rank
  colNames <- NULL
  switch(
    tax_rank,
    Strain = {
      colNames <- c("Genus", "Species", "Strain")
    },
    Species = {
      colNames <- c("Genus", "Species")
    },
    Genus = {
      colNames <- c("Genus")
    },
    Family = {
      colNames <- c("Family")
    },
    Order = {
      colNames <- c("Order")
    }
  )

  # Validate rank and return the corresponding columns
  if (!is.null(colNames)) {
    return(colNames)
  } else {
    stop("Please choose a valid taxonomy rank!", call. = FALSE)
  }
}


#' Parse Greengenes-style taxonomy strings including a Strain rank
#'
#' Helper parse function for \code{phyloseq::import_biom()} that converts a split
#' Greengenes taxonomy character vector into a named vector of taxonomic ranks.
#' In addition to standard Greengenes ranks, this parser includes an extra
#' \code{"Strain"} rank.
#'
#' The input \code{char.vec} is assumed to already be split into rank entries
#' (as \code{import_biom()} does for Greengenes-style taxonomy). Each entry is
#' expected to have a 3-character prefix such as \code{"k__"}, \code{"p__"}, etc.,
#' which is removed.
#'
#' @param char.vec Character vector of taxonomy rank strings (already split by rank),
#'   typically provided by \code{phyloseq::import_biom()}.
#'
#' @return A named character vector with names:
#'   \code{Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain}.
#'
#' @examples
#' \dontrun{
#' x <- c("k__Bacteria","p__Firmicutes","c__Bacilli","o__Lactobacillales",
#'        "f__Streptococcaceae","g__Streptococcus","s__aureus","t__USA300")
#' parse_taxonomy_strain(x)
#' }
#'
#' @keywords internal
parse_taxonomy_strain <- function(char.vec) {
  # Remove the Greengenes rank prefix (e.g., "k__", "p__") from each entry
  named.char.vec <- substring(char.vec, first = 4)

  # Assign taxonomic rank names (Greengenes ranks + Strain)
  names(named.char.vec) <- c("Kingdom", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species", "Strain")

  return(named.char.vec)
}

#' Extract an agglomerated feature table from a phyloseq BIOM object
#'
#' Agglomerates a \code{phyloseq} object at a requested taxonomic rank, then
#' constructs a feature table that combines:
#' \itemize{
#'   \item a single \code{taxonomy} label column (built by uniting selected taxonomy
#'     ranks), and
#'   \item the OTU/ASV abundance table.
#' }
#'
#' @param biom_object A \code{phyloseq} object (e.g., from \code{phyloseq::import_biom()}).
#' @param tax_rank Character. Taxonomic rank passed to \code{phyloseq::tax_glom()}
#'   (e.g., \code{"Genus"}, \code{"Species"}).
#' @param col_names Character vector of taxonomy columns to unite into a single
#'   \code{taxonomy} label (e.g., from \code{get_colNames_per_rank()}).
#'
#' @details
#' Agglomeration is performed with \code{NArm = TRUE}, which removes taxa with
#' missing assignments at the chosen \code{tax_rank}.
#'
#' The taxonomy label is built with \code{tidyr::unite(..., sep = "_")}, producing
#' labels such as \code{"Corynebacterium_propinq"} depending on your taxonomy table.
#'
#' @return A data.frame where the first column is \code{taxonomy} and remaining
#'   columns are sample abundances.
#'
#' @examples
#' \dontrun{
#' unite_cols <- get_colNames_per_rank("Species")
#' ft <- extract_table(biom_object, tax_rank = "Species", col_names = unite_cols)
#' head(ft)
#' }
#'
#' @export
extract_table <- function(biom_object, tax_rank, col_names) {
  # Agglomerate taxa at the requested taxonomic rank
  biom_object <- phyloseq::tax_glom(biom_object, taxrank = tax_rank, NArm = TRUE)

  # Build a single taxonomy label column and bind to the OTU/ASV table
  feature_table <- cbind(
    dplyr::select(
      tidyr::unite(
        data.frame(phyloseq::tax_table(biom_object)),
        taxonomy,
        all_of(col_names),
        sep = "_"
      ),
      "taxonomy"
    ),
    data.frame(phyloseq::otu_table(biom_object))
  )

  return(feature_table)
}

#' Clean a feature table by enforcing unique taxonomy labels and optional ordering
#'
#' Takes a feature table containing a \code{taxonomy} column and sample abundance
#' columns, makes taxonomy labels unique, converts the taxonomy column into rownames,
#' and optionally orders rows by decreasing mean abundance.
#'
#' @param feature_table Data frame with a \code{taxonomy} column and sample columns.
#' @param order_table Logical. If \code{TRUE}, order rows by decreasing \code{rowMeans()}.
#'
#' @details
#' \itemize{
#'   \item \code{make.unique(..., sep = "_")} ensures no duplicated taxonomy labels,
#'     which would otherwise break rowname-based indexing.
#'   \item Rownames are set using \code{tibble::column_to_rownames()}.
#'   \item If \code{order_table = TRUE}, rows are sorted by \code{rowMeans()} in
#'     decreasing order.
#' }
#'
#' @return A data frame with taxonomy in rownames and samples in columns. If
#'   \code{order_table = TRUE}, the rows are sorted by decreasing mean abundance.
#'
#' @examples
#' \dontrun{
#' ft_clean <- clean_table(feature_table, order_table = TRUE)
#' head(ft_clean)
#' }
#'
#' @export
clean_table <- function(feature_table, order_table) {
  # Ensure taxonomy labels are unique
  feature_table["taxonomy"] <- make.unique(feature_table$taxonomy, sep = "_")

  # Move taxonomy column into rownames
  feature_table <- tibble::column_to_rownames(
    tibble::remove_rownames(feature_table),
    var = "taxonomy"
  )

  # Optionally order features by mean abundance (high -> low)
  if (isTRUE(order_table)) {
    feature_table <- feature_table[order(rowMeans(feature_table), decreasing = TRUE), , drop = FALSE]
  }

  return(feature_table)
}


#' Filter low-abundance taxa based on maximum relative abundance
#'
#' Removes taxa (e.g., species/ASVs) whose relative abundance never reaches a
#' specified threshold in any sample. This is a simple prevalence/abundance
#' filter commonly used to reduce noise before ordination, clustering, or
#' differential abundance analyses.
#'
#' @param rel_abundance Numeric matrix or data.frame with taxa as rows and samples
#'   as columns, containing relative abundances (values typically in \code{[0, 1]}).
#' @param threshold Numeric. Minimum relative abundance required in at least one
#'   sample for a taxon to be retained (default \code{0.01}, i.e., 1%).
#'
#' @details
#' For each taxon (row), the maximum relative abundance across all samples is
#' computed. Taxa with \code{max(abundance) < threshold} are removed.
#'
#' Note that this filter:
#' \itemize{
#'   \item is sensitive to single-sample spikes (a taxon is kept if it exceeds the
#'     threshold in just one sample),
#'   \item does not account for prevalence (number of samples above threshold).
#' }
#' Consider combining with prevalence-based filters if needed.
#'
#' @return A matrix/data.frame of the same type as \code{rel_abundance}, containing
#'   only taxa whose maximum relative abundance meets or exceeds \code{threshold}.
#'
#' @examples
#' \dontrun{
#' rel <- matrix(
#'   c(0.2, 0.1, 0.0,
#'     0.005, 0.002, 0.001,
#'     0.03, 0.04, 0.01),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("Sp1","Sp2","Sp3"), c("S1","S2","S3"))
#' )
#'
#' filtered <- filter_low_abundance(rel, threshold = 0.01)
#' }
#'
#' @export
filter_low_abundance <- function(rel_abundance, threshold = 0.01) {
  # Compute maximum relative abundance per taxon across samples
  species_max <- apply(rel_abundance, 1, max)

  # Retain only taxa with max abundance >= threshold
  filtered_df <- rel_abundance[species_max >= threshold, , drop = FALSE]

  message(
    "Filtered from ", nrow(rel_abundance),
    " to ", nrow(filtered_df), " taxa"
  )

  return(filtered_df)
}


# ----- Within-timepoint replicate distances ----
# For replicate similarity
#' Align abundance and metadata, map timepoints, and normalize abundances per sample
#'
#' Prepares inputs for distance-based analyses by:
#' \itemize{
#'   \item intersecting samples between an abundance table and a metadata table,
#'   \item creating standardized metadata columns used downstream (sample/time/replicate/SynCom),
#'   \item mapping time labels (T1, T2, T3, TF) to numeric order (1..4), and
#'   \item producing a samples x taxa matrix normalized to relative abundances per sample.
#' }
#'
#' @param abund Numeric matrix/data.frame of taxa x samples (columns are sample IDs).
#' @param meta Data frame of sample metadata with rownames as sample IDs and at least
#'   \code{ATTRIBUTE_SynCom}, \code{ATTRIBUTE_Time}, and \code{ATTRIBUTE_Replicate}.
#'
#' @details
#' **Sample alignment**
#' Samples are kept in the intersection of \code{colnames(abund)} and \code{rownames(meta)}.
#' Metadata is reordered to match the abundance columns.
#'
#' **Time mapping**
#' \code{ATTRIBUTE_Time} must be one of \code{T1}, \code{T2}, \code{T3}, \code{TF}.
#' These are mapped to numeric order \code{1..4} and also stored as an ordered factor.
#'
#' **Normalization**
#' The abundance matrix is transposed to samples x taxa and normalized by row sums
#' (per-sample relative abundances). NAs are set to 0 and all-zero rows are protected
#' by replacing row sums of 0 with 1 before division.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{meta}: aligned and augmented metadata data.frame.
#'   \item \code{X}: numeric matrix of samples x taxa (relative abundances).
#' }
#'
#' @examples
#' \dontrun{
#' prep <- prepare_data(abund = rel_abund_table, meta = md)
#' head(prep$meta)
#' dim(prep$X)  # samples x taxa
#' }
#'
#' @export
prepare_data <- function(abund, meta) {
  # Ensure column/rownames are present
  stopifnot(!is.null(colnames(abund)), !is.null(rownames(meta)))

  # Intersect samples
  common_samples <- intersect(colnames(abund), rownames(meta))
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between abund colnames and meta rownames.")
  }

  abund2 <- abund[, common_samples, drop = FALSE]
  meta2  <- meta[common_samples, , drop = FALSE]

  # Standardize key fields used downstream
  meta2 <- meta2 %>%
    mutate(
      sample_id = rownames(meta2),
      syncom_id = .data$ATTRIBUTE_SynCom,
      time_raw  = .data$ATTRIBUTE_Time,
      rep_raw   = .data$ATTRIBUTE_Replicate
    )

  # Map times to numeric order 1..4 (T1,T2,T3,TF->4)
  time_map <- c(T1 = 1, T2 = 2, T3 = 3, TF = 4)
  if (!all(meta2$time_raw %in% names(time_map))) {
    bad_levels <- setdiff(unique(meta2$time_raw), names(time_map))
    stop("Unexpected time labels in ATTRIBUTE_Time: ", paste(bad_levels, collapse = ", "))
  }

  meta2 <- meta2 %>%
    mutate(
      time_num   = unname(time_map[time_raw]),
      time_label = factor(time_raw, levels = c("T1","T2","T3","TF"))
    )

  # Normalize abundances per sample (samples x taxa)
  X <- t(abund2)
  X[is.na(X)] <- 0
  row_sums <- rowSums(X)
  row_sums[row_sums == 0] <- 1
  X <- sweep(X, 1, row_sums, "/")

  list(meta = meta2, X = as.matrix(X))
}

#' Compute within-timepoint replicate dissimilarities per SynCom
#'
#' For each SynCom and each time point, computes all pairwise dissimilarities among
#' replicate samples (e.g., R1 vs R2, R1 vs R3, R2 vs R3) using \code{vegan::vegdist()}.
#' Returns a tidy table containing individual pairwise distances plus per-group
#' summary statistics (mean, SD, number of pairs).
#'
#' @param meta Metadata data.frame produced by \code{prepare_data()}, containing at least:
#'   \code{sample_id}, \code{syncom_id}, \code{time_num}, \code{time_label}, \code{rep_raw}.
#' @param X Numeric matrix of samples x taxa (rows are sample IDs), typically \code{prepare_data()$X}.
#' @param method Character. Distance method passed to \code{vegan::vegdist()}
#'   (default \code{"bray"}).
#'
#' @details
#' Distances are computed within blocks defined by \code{syncom_id x time_num}.
#' If a block has fewer than 2 samples, no rows are returned for that block.
#'
#' The output includes:
#' \itemize{
#'   \item \code{replicate_pair}: label like \code{"R1_vs_R2"} derived from \code{rep_raw}
#'   \item \code{pairwise_dist}: the computed dissimilarity
#'   \item \code{mean_dist}, \code{sd_dist}, \code{n_pairs}: repeated per block for convenience
#' }
#'
#' @return A tibble with one row per replicate pair within each SynCom-timepoint block,
#'   plus summary columns per block.
#'
#' @examples
#' \dontrun{
#' prep <- prepare_data(abund, meta)
#' dist_tbl <- compute_within_tp_distances(prep$meta, prep$X, method = "bray")
#' dist_tbl
#' }
#'
#' @export
compute_within_tp_distances <- function(meta, X, method = "bray") {
  dat <- meta %>% select(sample_id, syncom_id, time_num, time_label, rep_raw)

  # Helper to compute pairwise distances for a block of sample_ids
  compute_block <- function(sample_ids) {
    if (length(sample_ids) < 2) return(tibble())

    # Distances among samples (rows)
    d  <- vegdist(X[sample_ids, , drop = FALSE], method = method)
    dv <- as.numeric(d)

    # Label replicate pairs (e.g., "R1_vs_R2")
    reps <- meta$rep_raw[match(sample_ids, meta$sample_id)]
    pair_names <- combn(reps, 2, FUN = function(xx) paste0(xx[1], "_vs_", xx[2]))

    tibble(replicate_pair = pair_names, pairwise_dist = dv)
  }

  # Group by SynCom x Time and compute the pairwise distances per block
  dist_tbl <- dat %>%
    group_by(syncom_id, time_num, time_label) %>%
    group_modify(~{
      sample_ids <- .x$sample_id
      compute_block(sample_ids)
    }) %>%
    ungroup()

  # Add per-block summary statistics while retaining individual pairs
  dist_tbl %>%
    group_by(syncom_id, time_num, time_label) %>%
    mutate(
      mean_dist = mean(pairwise_dist),
      sd_dist   = sd(pairwise_dist),
      n_pairs   = dplyr::n()
    ) %>%
    ungroup()
}

#' Plot within-timepoint replicate dissimilarity trajectories per SynCom
#'
#' Visualizes replicate dissimilarities over time with one facet per SynCom.
#' Individual replicate-pair distances are shown as jittered points; the mean per
#' time point is overlaid as a line and points (via \code{stat_summary}).
#'
#' @param dist_tbl Tibble produced by \code{compute_within_tp_distances()}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_replicate_similarity(dist_tbl)
#' p
#' }
#'
#' @export
plot_replicate_similarity <- function(dist_tbl) {
  # Order SynCom facets numerically based on digits in syncom_id (e.g., "SynCom10")
  dist_tbl <- dist_tbl %>%
    mutate(
      syncom_order = as.numeric(gsub("\\D", "", syncom_id)),
      syncom_id = factor(syncom_id, levels = unique(syncom_id[order(syncom_order)]))
    )

  ggplot(dist_tbl, aes(x = time_num, y = pairwise_dist)) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.05, height = 0)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    facet_wrap(~ syncom_id, scales = "free_y") +
    scale_x_continuous(breaks = 1:4, labels = c("T1","T2","T3","TF")) +
    labs(
      x = "Time point",
      y = "Replicate dissimilarity (Bray-Curtis)",
      title = "Within-time point replicate dissimilarity per SynCom"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


#' Compute distance to the final time point (TF) within each SynCom
#'
#' Quantifies "stabilization" by measuring, for each sample, its dissimilarity to
#' the SynCom-specific final state (TF). Distances are computed within each SynCom
#' using \code{vegan::vegdist()} (default Bray-Curtis), and summarized per SynCom
#' and time point.
#'
#' Three modes are supported:
#' \itemize{
#'   \item \code{"centroid"}: distance from each sample to the TF centroid (mean profile).
#'   \item \code{"replicate"}: distance from each sample to the TF sample(s) with the
#'     same replicate ID; falls back to centroid if no match exists.
#'   \item \code{"allpair_mean"}: mean distance from each sample to all TF replicates.
#' }
#'
#' @param meta Metadata data frame (typically \code{prepare_data()$meta}) with sample IDs
#'   and SynCom/time/replicate information. Must correspond row-for-row to \code{X}.
#' @param X Numeric matrix of samples x taxa (rows are samples). Must correspond row-for-row
#'   to \code{meta}. Typically \code{prepare_data()$X}.
#' @param method Character. Distance method passed to \code{vegan::vegdist()}
#'   (default \code{"bray"}).
#' @param mode Character. How to define the TF reference:
#'   \code{"centroid"}, \code{"replicate"}, or \code{"allpair_mean"}.
#'
#' @details
#' **Required metadata fields**
#' The function expects \code{meta} to include:
#' \code{sample_id}, \code{syncom_id}, \code{time_label}, \code{time_num}, \code{rep_raw}.
#' If any are missing, it attempts to create them from:
#' \code{ATTRIBUTE_SynCom}, \code{ATTRIBUTE_Time}, \code{ATTRIBUTE_Replicate}.
#'
#' **TF detection**
#' TF samples are those with \code{time_label == "TF"} within each SynCom. If a SynCom
#' has no TF samples, distances for that SynCom are returned as \code{NA}.
#'
#' **Computation notes**
#' \itemize{
#'   \item \code{"centroid"} uses \code{colMeans()} across TF samples as the reference vector.
#'   \item \code{"replicate"} matches TF sample(s) by replicate ID; if multiple matches are
#'     found, their distances are averaged; if none, the TF centroid is used.
#'   \item \code{"allpair_mean"} averages the distances to all TF replicates.
#' }
#'
#' @return A list with:
#' \itemize{
#'   \item \code{per_sample}: tibble with one row per sample, including \code{dist_to_TF}
#'     and the key metadata columns.
#'   \item \code{summary}: tibble summarized per \code{syncom_id x time_label x time_num}
#'     with mean, SD, and sample counts.
#' }
#'
#' @examples
#' \dontrun{
#' prep <- prepare_data(abund, meta)
#' out <- compute_distance_to_final(prep$meta, prep$X, method = "bray", mode = "centroid")
#' head(out$per_sample)
#' out$summary
#' }
#'
#' @export
compute_distance_to_final <- function(meta, X, method = "bray",
                                      mode = c("centroid", "replicate", "allpair_mean")) {
  mode <- match.arg(mode)
  stopifnot(nrow(meta) == nrow(X))

  # Coerce to tibble to avoid rowname weirdness
  meta <- tibble::as_tibble(meta)

  # ---- Ensure required columns exist (or derive them from ATTRIBUTE_* columns) ----
  if (!"sample_id" %in% names(meta)) {
    stop("meta must contain a 'sample_id' column. Did you pass prepared$meta?")
  }

  if (!"syncom_id" %in% names(meta)) {
    if ("ATTRIBUTE_SynCom" %in% names(meta)) {
      meta <- meta %>% mutate(syncom_id = .data$ATTRIBUTE_SynCom)
    } else {
      stop("meta lacks 'syncom_id' and 'ATTRIBUTE_SynCom'.")
    }
  }

  if (!"time_label" %in% names(meta) || !"time_num" %in% names(meta)) {
    if ("ATTRIBUTE_Time" %in% names(meta)) {
      time_map <- c(T1 = 1, T2 = 2, T3 = 3, TF = 4)
      meta <- meta %>%
        mutate(
          time_label = factor(.data$ATTRIBUTE_Time, levels = c("T1","T2","T3","TF")),
          time_num   = unname(time_map[as.character(.data$ATTRIBUTE_Time)])
        )
    } else {
      stop("meta lacks 'time_label/time_num' and 'ATTRIBUTE_Time'.")
    }
  }

  if (!"rep_raw" %in% names(meta)) {
    if ("ATTRIBUTE_Replicate" %in% names(meta)) {
      meta <- meta %>% mutate(rep_raw = .data$ATTRIBUTE_Replicate)
    } else {
      stop("meta lacks 'rep_raw' and 'ATTRIBUTE_Replicate'.")
    }
  }

  # Compact metadata used downstream
  dat <- meta %>%
    dplyr::select(sample_id, syncom_id, time_label, time_num, rep_raw)

  # NOTE: This helper is currently unused in the implementation, but kept as a
  # convenience for future refactoring.
  bc_to_vec <- function(sample_row, ref_vec) {
    d <- vegan::vegdist(rbind(sample_row, ref_vec), method = method)
    as.numeric(d)[1]
  }

  # Debug preview of the metadata block structure
  print(head(dat))

  # ---- Compute distance-to-TF within each SynCom ----
  results <- dat %>%
    dplyr::group_by(syncom_id) %>%
    dplyr::group_modify(~{
      block_meta <- .x
      ids   <- block_meta$sample_id
      tf_ids <- block_meta$sample_id[block_meta$time_label == "TF"]

      # If no TF samples exist for this SynCom, return NA distances
      if (length(tf_ids) == 0) {
        return(tibble::tibble(sample_id = ids, dist_to_TF = NA_real_))
      }

      if (mode == "centroid") {
        # Reference = TF centroid (mean profile)
        tf_centroid <- colMeans(X[tf_ids, , drop = FALSE])

        tibble::tibble(
          sample_id = ids,
          dist_to_TF = apply(X[ids, , drop = FALSE], 1, function(r) {
            d <- vegan::vegdist(rbind(r, tf_centroid), method = method)
            as.numeric(d)[1]
          })
        )

      } else if (mode == "replicate") {
        # Reference = TF sample(s) matching replicate ID; fallback to centroid
        tf_centroid <- colMeans(X[tf_ids, , drop = FALSE])

        tf_rep_map <- tibble::tibble(
          tf_id  = tf_ids,
          tf_rep = block_meta$rep_raw[match(tf_ids, block_meta$sample_id)]
        )

        purrr::map_dfr(ids, function(sid) {
          rep_s <- block_meta$rep_raw[block_meta$sample_id == sid]
          tf_match <- tf_rep_map$tf_id[tf_rep_map$tf_rep == rep_s]

          dval <- if (length(tf_match) >= 1) {
            if (length(tf_match) == 1) {
              as.numeric(vegan::vegdist(rbind(X[sid, ], X[tf_match, ]), method = method))[1]
            } else {
              # If multiple TF matches for the replicate, average distances
              mean(as.numeric(vegan::vegdist(rbind(X[sid, ], X[tf_match, ]), method = method))[seq_along(tf_match)])
            }
          } else {
            # Fallback: distance to TF centroid
            as.numeric(vegan::vegdist(rbind(X[sid, ], tf_centroid), method = method))[1]
          }

          tibble::tibble(sample_id = sid, dist_to_TF = dval)
        })

      } else { # mode == "allpair_mean"
        # Reference = average distance to all TF replicates
        purrr::map_dfr(ids, function(sid) {
          d <- vegan::vegdist(rbind(X[sid, ], X[tf_ids, ]), method = method)
          tibble::tibble(
            sample_id = sid,
            dist_to_TF = mean(as.numeric(d)[seq_along(tf_ids)])
          )
        })
      }
    }) %>%
    dplyr::ungroup() %>%
    # Drop the group key that group_modify appended and reattach clean metadata
    dplyr::select(sample_id, dist_to_TF) %>%
    dplyr::left_join(dat, by = "sample_id")

  # ---- Summarize per SynCom x time point ----
  summary_tbl <- results %>%
    dplyr::group_by(syncom_id, time_label, time_num) %>%
    dplyr::summarize(
      mean_dist_to_TF = mean(dist_to_TF, na.rm = TRUE),
      sd_dist_to_TF   = sd(dist_to_TF, na.rm = TRUE),
      n               = sum(!is.na(dist_to_TF)),
      .groups = "drop"
    )

  list(per_sample = results, summary = summary_tbl)
}

#' Plot distance-to-final-state trajectories per SynCom
#'
#' Creates a faceted plot (one panel per SynCom) showing per-sample Bray-Curtis
#' distance to the SynCom-specific final state (TF) over time. Individual samples
#' are shown as jittered points, and the mean distance at each time point is
#' overlaid as a line and points.
#'
#' @param per_sample Tibble produced by \code{compute_distance_to_final()$per_sample}.
#'   Must include \code{syncom_id}, \code{time_num}, and \code{dist_to_TF}.
#' @param summary_tbl Tibble produced by \code{compute_distance_to_final()$summary}.
#'   Must include \code{syncom_id}, \code{time_num}, and \code{mean_dist_to_TF}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' out <- compute_distance_to_final(meta, X, mode = "allpair_mean")
#' p <- plot_distance_to_final(out$per_sample, out$summary)
#' p
#' }
#'
#' @export
plot_distance_to_final <- function(per_sample, summary_tbl) {
  # Consistent numeric ordering of SynCom facet levels (based on digits in syncom_id)
  order_levels <- summary_tbl %>%
    mutate(syncom_order = as.numeric(gsub("\\D", "", syncom_id))) %>%
    arrange(syncom_order) %>%
    pull(syncom_id) %>%
    unique()

  per_sample <- per_sample %>%
    mutate(syncom_id = factor(syncom_id, levels = order_levels))
  summary_tbl <- summary_tbl %>%
    mutate(syncom_id = factor(syncom_id, levels = order_levels))

  ggplot() +
    geom_point(
      data = per_sample,
      aes(x = time_num, y = dist_to_TF),
      alpha = 0.6,
      position = position_jitter(width = 0.05, height = 0)
    ) +
    geom_line(
      data = summary_tbl,
      aes(x = time_num, y = mean_dist_to_TF, group = 1),
      linewidth = 0.8
    ) +
    geom_point(
      data = summary_tbl,
      aes(x = time_num, y = mean_dist_to_TF),
      size = 2
    ) +
    facet_wrap(~ syncom_id, scales = "free_y") +
    scale_x_continuous(breaks = 1:4, labels = c("T1","T2","T3","TF")) +
    labs(
      x = "Time point",
      y = "Bray-Curtis distance to final state (TF)",
      title = "Stabilization toward final community composition"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}
