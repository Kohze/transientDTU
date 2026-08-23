#' Legacy list representation of the synthetic example
#'
#' A known-truth simulation with 12 two-feature genes, three groups, eight
#' ordered stages, and two replicates per group-stage cell. Genes represent
#' null, monotonic, persistent-divergence, and bounded transient trajectories.
#' It is small enough for examples and tests and contains no biological or
#' personally identifiable data. New code should generally use
#' [transientExampleSE] and [transientExamplePairwise], which follow standard
#' Bioconductor container conventions.
#'
#' @format A named list with seven principal elements:
#' \describe{
#'   \item{pairwise}{A `data.frame` of directed focal-versus-comparator effects,
#'     raw p-values, globally adjusted p-values, and gene-level adjusted
#'     p-values.}
#'   \item{usage}{A 24-feature by 48-sample proportion matrix.}
#'   \item{col_data}{Sample-level group, stage, and replicate metadata.}
#'   \item{truth}{The true feature-level bounded episodes.}
#'   \item{gene_truth}{The trajectory class assigned to every gene.}
#'   \item{stage_order}{The complete ordered stage labels.}
#'   \item{groups}{The three group labels.}
#' }
#'
#' @source Generated entirely by [simulateTransientDTU()] with seed 20260823;
#'   see `inst/scripts/make-example-data.R`.
#'
#' @usage data("transientExample")
#' @examples
#' data("transientExample")
#' table(transientExample$gene_truth$trajectory)
#' dim(transientExample$usage)
#'
"transientExample"

#' Synthetic transcript usage in a SummarizedExperiment
#'
#' The primary Bioconductor representation of the package's reproducible
#' known-truth example. Rows are transcript features, columns are samples, and
#' the `usage` assay contains transcript-usage proportions. Feature annotation
#' is stored in `rowData`, sample design variables are stored in `colData`, and
#' the declared stage order, group labels, simulation truth, parameters, and
#' provenance are stored in `metadata`.
#'
#' @format A [SummarizedExperiment::SummarizedExperiment] with 24 rows and 48
#'   columns:
#' \describe{
#'   \item{assay `usage`}{Feature-by-sample transcript-usage proportions.}
#'   \item{`rowData`}{`feature_id`, `gene_id`, and `gene_name`.}
#'   \item{`colData`}{`sample_id`, `group`, `stage`, `replicate`, and
#'     `stage_index`.}
#'   \item{`metadata`}{`stage_order`, `groups`, feature- and gene-level truth,
#'     simulation parameters, and a provenance statement.}
#' }
#'
#' @source Generated entirely by [simulateTransientDTU()] with seed 20260823;
#'   see `inst/scripts/make-example-data.R`.
#'
#' @usage data("transientExampleSE")
#' @examples
#' data("transientExampleSE")
#' SummarizedExperiment::assayNames(transientExampleSE)
#' head(SummarizedExperiment::colData(transientExampleSE))
#' S4Vectors::metadata(transientExampleSE)$stage_order
#'
"transientExampleSE"

#' Synthetic directed pairwise DTU evidence
#'
#' Generic long-format, directed focal-versus-comparator evidence matched to
#' [transientExampleSE]. It represents the standard table contract accepted by
#' [auditTransientInput()], [makeStageDTU()], and [runTransientDTU()]. The
#' signed effect is focal usage minus comparator usage.
#'
#' @format An [S4Vectors::DataFrame] with one row per feature, focal group,
#'   comparator group, and stage. It contains `feature_id`, `gene_id`,
#'   `gene_name`, `focal_group`, `comparator_group`, `stage`, `effect`,
#'   `p_value`, `q_value`, and `gene_q`.
#'
#' @source Generated entirely by [simulateTransientDTU()] with seed 20260823;
#'   see `inst/scripts/make-example-data.R`.
#'
#' @usage data("transientExamplePairwise")
#' @examples
#' data("transientExamplePairwise")
#' head(transientExamplePairwise)
#'
"transientExamplePairwise"
