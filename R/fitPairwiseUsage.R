#' Fit stage-specific pairwise usage contrasts with limma
#'
#' Convenience inference adapter for users starting from replicate-level usage
#' proportions. It fits a saturated group-by-stage cell-means model with
#' [limma::lmFit()], evaluates every directed focal-versus-comparator contrast
#' within each stage, and returns the generic table accepted by
#' [makeStageDTU()]. The detector itself remains independent of limma.
#'
#' @param usage A feature-by-sample proportion matrix or a
#'   [SummarizedExperiment::SummarizedExperiment].
#' @param col_data Sample metadata for matrix input.
#' @param feature_data Feature metadata for matrix input. Its row names must
#'   match usage rows, unless `feature_col` identifies them explicitly.
#' @param assay_name Assay used for `SummarizedExperiment` input.
#' @param sample_col Optional sample identifier in `col_data`.
#' @param feature_col Optional feature identifier in `feature_data`.
#' @param gene_col Gene identifier in feature metadata.
#' @param gene_name_col Optional gene display-name column.
#' @param group_col,stage_col Sample metadata columns.
#' @param stage_order Complete ordered stages. Inferred from a factor's levels
#'   or first appearance when omitted; supplying it explicitly is recommended.
#' @param adjustment Adjustment method applied across the complete directed
#'   feature-by-stage-by-pair family.
#' @param trend,robust Arguments passed to [limma::eBayes()].
#'
#' @return A base `data.frame` satisfying the default [makeStageDTU()] schema.
#'
#' @details
#' Gene-level p-values are Simes combinations of all raw component p-values for
#' a gene and are then adjusted across genes. They are provided as an optional
#' screen, not as episode-level error control. Direct proportion-scale limma is
#' a convenience adapter; users with a preferred DTU model should supply its
#' results directly to [makeStageDTU()].
#'
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'     simulated <- simulateTransientDTU(n_genes = 8, seed = 12)
#'     pairwise <- fitPairwiseUsage(
#'         simulated$usage,
#'         simulated$col_data,
#'         feature_data = data.frame(
#'             gene_id = sub("_tx[12]$", "", rownames(simulated$usage)),
#'             row.names = rownames(simulated$usage)
#'         ),
#'         stage_order = simulated$stage_order
#'     )
#'     head(pairwise)
#' }
#'
#' @export
fitPairwiseUsage <- function(
    usage,
    col_data = NULL,
    feature_data = NULL,
    assay_name = "usage",
    sample_col = NULL,
    feature_col = NULL,
    gene_col = "gene_id",
    gene_name_col = NULL,
    group_col = "group",
    stage_col = "stage",
    stage_order = NULL,
    adjustment = "BH",
    trend = TRUE,
    robust = TRUE
) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Install the Bioconductor package 'limma' to use this adapter.",
            call. = FALSE)
    }
    if (!adjustment %in% stats::p.adjust.methods) {
        stop("Unknown multiplicity adjustment method.", call. = FALSE)
    }
    if (methods::is(usage, "SummarizedExperiment")) {
        extracted <- .extract_usage(
            usage, assay_name = assay_name,
            group_col = group_col, stage_col = stage_col
        )
        feature_metadata <- as.data.frame(
            SummarizedExperiment::rowData(usage)
        )
        rownames(feature_metadata) <- rownames(usage)
    } else {
        extracted <- .extract_usage(
            usage, col_data = col_data, sample_col = sample_col,
            group_col = group_col, stage_col = stage_col
        )
        if (is.null(feature_data)) {
            stop("'feature_data' is required when 'usage' is a matrix.",
                call. = FALSE)
        }
        feature_metadata <- .assert_data_frame(feature_data, "feature_data")
    }
    matrix_value <- extracted$usage
    metadata <- extracted$col_data
    .assert_columns(feature_metadata, gene_col, "feature_data")
    feature_ids <- if (is.null(feature_col)) {
        rownames(feature_metadata)
    } else {
        .assert_columns(feature_metadata, feature_col, "feature_data")
        as.character(feature_metadata[[feature_col]])
    }
    if (is.null(feature_ids) || anyNA(feature_ids) ||
            anyDuplicated(feature_ids)) {
        stop(
            "Feature metadata needs unique feature identifiers.",
            call. = FALSE
        )
    }
    matched_features <- match(rownames(matrix_value), feature_ids)
    if (anyNA(matched_features)) {
        stop("Every usage row must match one feature-metadata row.",
            call. = FALSE)
    }
    feature_metadata <- feature_metadata[matched_features, , drop = FALSE]
    feature_ids <- rownames(matrix_value)
    gene_ids <- .as_id(feature_metadata[[gene_col]], gene_col)
    gene_names <- if (is.null(gene_name_col)) {
        gene_ids
    } else {
        .assert_columns(feature_metadata, gene_name_col, "feature_data")
        .as_id(
            feature_metadata[[gene_name_col]], gene_name_col, allow_na = TRUE
        )
    }

    groups <- sort(unique(as.character(metadata[[group_col]])))
    observed_stages <- as.character(metadata[[stage_col]])
    if (is.null(stage_order)) {
        stage_value <- metadata[[stage_col]]
        stage_order <- if (is.factor(stage_value)) {
            levels(stage_value)
        } else {
            unique(observed_stages)
        }
    }
    stage_order <- as.character(stage_order)
    if (any(!observed_stages %in% stage_order) || anyDuplicated(stage_order)) {
        stop("'stage_order' must uniquely include every observed stage.",
            call. = FALSE)
    }
    cell <- interaction(
        factor(metadata[[group_col]], levels = groups),
        factor(metadata[[stage_col]], levels = stage_order),
        drop = TRUE, sep = "__"
    )
    expected_cells <- as.vector(outer(groups, stage_order, paste, sep = "__"))
    missing_cells <- setdiff(expected_cells, levels(cell))
    if (length(missing_cells)) {
        stop(
            "The saturated design lacks group-stage cells: ",
            paste(missing_cells, collapse = ", "), ".",
            call. = FALSE
        )
    }
    cell <- factor(as.character(cell), levels = expected_cells)
    design <- stats::model.matrix(~ 0 + cell)
    colnames(design) <- expected_cells
    if (qr(design)$rank != ncol(design)) {
        stop("The group-stage design matrix is not full rank.", call. = FALSE)
    }

    pair_metadata <- do.call(rbind, lapply(stage_order, function(stage) {
        do.call(rbind, lapply(groups, function(focal) {
            data.frame(
                stage = stage,
                focal_group = focal,
                comparator_group = setdiff(groups, focal),
                stringsAsFactors = FALSE
            )
        }))
    }))
    contrast <- matrix(
        0,
        nrow = ncol(design), ncol = nrow(pair_metadata),
        dimnames = list(
            colnames(design),
            paste(
                pair_metadata$focal_group,
                pair_metadata$comparator_group,
                pair_metadata$stage,
                sep = "__vs__"
            )
        )
    )
    for (index in seq_len(nrow(pair_metadata))) {
        focal_cell <- paste(
            pair_metadata$focal_group[[index]],
            pair_metadata$stage[[index]], sep = "__"
        )
        comparator_cell <- paste(
            pair_metadata$comparator_group[[index]],
            pair_metadata$stage[[index]], sep = "__"
        )
        contrast[focal_cell, index] <- 1
        contrast[comparator_cell, index] <- -1
    }
    fit <- limma::lmFit(matrix_value, design)
    fit <- limma::contrasts.fit(fit, contrast)
    fit <- limma::eBayes(fit, trend = trend, robust = robust)
    effects <- fit$coefficients
    p_values <- fit$p.value
    q_values <- matrix(
        stats::p.adjust(as.vector(p_values), method = adjustment),
        nrow = nrow(p_values), ncol = ncol(p_values),
        dimnames = dimnames(p_values)
    )

    gene_p <- vapply(split(seq_along(gene_ids), gene_ids), function(index) {
        values <- sort(as.vector(p_values[index, , drop = FALSE]))
        min(1, min(values * length(values) / seq_along(values)))
    }, numeric(1))
    gene_q <- stats::p.adjust(gene_p, method = adjustment)

    rows <- lapply(seq_len(ncol(effects)), function(index) {
        data.frame(
            feature_id = feature_ids,
            gene_id = gene_ids,
            gene_name = gene_names,
            focal_group = pair_metadata$focal_group[[index]],
            comparator_group = pair_metadata$comparator_group[[index]],
            stage = pair_metadata$stage[[index]],
            effect = effects[, index],
            p_value = p_values[, index],
            q_value = q_values[, index],
            gene_q = unname(gene_q[gene_ids]),
            stringsAsFactors = FALSE
        )
    })
    answer <- do.call(rbind, rows)
    rownames(answer) <- NULL
    answer
}
