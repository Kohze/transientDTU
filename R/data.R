#' Synthetic ordered multi-group transcript-usage example
#'
#' A known-truth simulation with 12 two-feature genes, three groups, eight
#' ordered stages, and two replicates per group-stage cell. Genes represent
#' null, monotonic, persistent-divergence, and bounded transient trajectories.
#' It is small enough for examples and vignettes and contains no biological or
#' personally identifiable data.
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
#'   see `data-raw/make-example-data.R` in the source repository.
#'
#' @usage data("transientExample")
#' @examples
#' data("transientExample")
#' table(transientExample$gene_truth$trajectory)
#' dim(transientExample$usage)
#'
"transientExample"
