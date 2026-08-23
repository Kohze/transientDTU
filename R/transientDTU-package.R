#' transientDTU: bounded temporal transcript-usage patterns
#'
#' `transientDTU` is a post-inference framework for detecting explicit
#' diverge-reconverge episodes in ordered, multi-group transcript-usage
#' experiments. It consumes generic pairwise DTU evidence rather than replacing
#' upstream statistical models. The package formalizes stage eligibility,
#' episode boundaries, replicate separation, reciprocal exchanges, and stable
#' candidate ranking, while exposing simulation and sensitivity tools needed to
#' audit the decision layer.
#'
#' The primary workflow is:
#'
#' 1. [makeStageDTU()] standardizes pairwise evidence;
#' 2. [detectEpisodes()] identifies bounded runs;
#' 3. [checkReplicateSeparation()] optionally checks non-overlap;
#' 4. [annotateReciprocal()] marks compositional exchanges;
#' 5. [rankCandidates()] creates a deterministic panel; or
#' 6. [runTransientDTU()] performs all steps together.
#'
#' [SummarizedExperiment::SummarizedExperiment] assays can be used for
#' replicate checks and the optional [fitPairwiseUsage()] inference adapter.
#'
#' @references
#' Van den Berge K, Soneson C, Robinson MD, Clement L (2017). stageR: a general
#' stage-wise method for controlling the gene-level false discovery rate in
#' differential expression and differential transcript usage. *Genome Biology*,
#' 18, 151. \doi{10.1186/s13059-017-1277-0}.
#'
#' Gilis J, Vitting-Seerup K, Van den Berge K, Clement L (2023). satuRn:
#' scalable analysis of differential transcript usage for bulk and single-cell
#' RNA-sequencing applications. *Genome Biology*, 24, 24.
#' \doi{10.1186/s13059-023-02863-7}.
#'
#' @keywords internal
"_PACKAGE"
