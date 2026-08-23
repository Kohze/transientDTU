.audit_issue <- function(severity, code, message, n_rows = NA_integer_) {
    data.frame(
        severity = severity,
        code = code,
        n_rows = as.integer(n_rows),
        message = message,
        stringsAsFactors = FALSE
    )
}

.empty_audit_coverage <- function() {
    data.frame(
        focal_group = character(),
        stage = character(),
        expected_feature_cells = integer(),
        observed_feature_cells = integer(),
        complete_feature_cells = integer(),
        completeness_fraction = numeric(),
        observed_atomic_rows = integer(),
        expected_atomic_rows = integer(),
        stringsAsFactors = FALSE
    )
}

.as_audit_numeric <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    value <- trimws(as.character(x))
    answer <- rep(NA_real_, length(value))
    numeric_pattern <- paste0(
        "^[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)",
        "(?:[eE][+-]?[0-9]+)?$"
    )
    valid <- !is.na(value) & grepl(numeric_pattern, value, perl = TRUE)
    answer[valid] <- as.numeric(value[valid])
    infinite <- match(value, c("Inf", "+Inf", "-Inf"), nomatch = 0L)
    answer[infinite == 1L | infinite == 2L] <- Inf
    answer[infinite == 3L] <- -Inf
    answer
}

.finish_input_audit <- function(summary, issues, coverage, replicate_cells) {
    issue_table <- if (length(issues)) {
        do.call(rbind, issues)
    } else {
        data.frame(
            severity = character(), code = character(),
            n_rows = integer(), message = character(),
            stringsAsFactors = FALSE
        )
    }
    status <- if (any(issue_table$severity == "error")) {
        "fail"
    } else if (any(issue_table$severity == "warning")) {
        "review"
    } else {
        "ready"
    }
    summary <- rbind(
        data.frame(
            metric = "status", value = status,
            stringsAsFactors = FALSE
        ),
        summary
    )
    structure(
        list(
            status = status,
            ready = identical(status, "ready"),
            summary = S4Vectors::DataFrame(summary),
            issues = S4Vectors::DataFrame(issue_table),
            coverage = S4Vectors::DataFrame(coverage),
            replicate_cells = S4Vectors::DataFrame(replicate_cells)
        ),
        class = c("TransientDTUAudit", "list")
    )
}

#' Audit a proposed transient-DTU input before analysis
#'
#' Performs structural and coverage checks without applying discovery
#' thresholds. The audit identifies blocking schema errors, incomplete directed
#' contrast families, missing inferential values, and replicate-data gaps. It
#' cannot verify effect orientation, the upstream model, or the multiplicity
#' family, so those remain explicit informational checks.
#'
#' @param pairwise Proposed generic pairwise evidence table.
#' @param stage_order Complete biological stage order.
#' @param feature_col,gene_col,focal_col,comparator_col Identifier columns.
#' @param stage_col,effect_col,q_col Evidence columns.
#' @param gene_name_col,gene_q_col Optional annotation and gene-level q-value
#'   columns.
#' @param expected_groups Complete group set; inferred when `NULL`.
#' @param expected_focal_groups Groups intended to act as focal groups. By
#'   default these are inferred from the focal column, allowing targeted scans.
#' @param min_comparators Minimum comparators required for the intended rule.
#'   Use one for a two-group design.
#' @param require_all_comparators Whether every non-focal expected group is
#'   required.
#' @param usage Optional usage matrix or
#'   [SummarizedExperiment::SummarizedExperiment] for replicate-data auditing.
#' @param col_data Sample metadata for matrix input.
#' @param assay_name,sample_col,group_col Metadata controls passed to the usage
#'   extractor.
#' @param sample_stage_col Stage column in sample metadata.
#'
#' @return A `TransientDTUAudit` list with `status`, `ready`, `summary`,
#'   `issues`, `coverage`, and `replicate_cells`. Tables are
#'   [S4Vectors::DataFrame] objects.
#'
#' @details
#' `status = "fail"` indicates at least one structural error; `"review"`
#' indicates non-blocking warnings; and `"ready"` means that only informational
#' limitations remain. A ready audit is not evidence that upstream inference or
#' multiplicity adjustment was scientifically appropriate.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 12, seed = 15)
#' audit <- auditTransientInput(
#'     simulated$pairwise,
#'     simulated$stage_order,
#'     usage = simulated$usage,
#'     col_data = simulated$col_data
#' )
#' audit
#' head(audit$coverage)
#'
#' @export
auditTransientInput <- function(
    pairwise,
    stage_order,
    feature_col = "feature_id",
    gene_col = "gene_id",
    focal_col = "focal_group",
    comparator_col = "comparator_group",
    stage_col = "stage",
    effect_col = "effect",
    q_col = "q_value",
    gene_name_col = NULL,
    gene_q_col = NULL,
    expected_groups = NULL,
    expected_focal_groups = NULL,
    min_comparators = 2L,
    require_all_comparators = TRUE,
    usage = NULL,
    col_data = NULL,
    assay_name = "usage",
    sample_col = NULL,
    group_col = "group",
    sample_stage_col = "stage"
) {
    pairwise <- .assert_data_frame(pairwise, "pairwise")
    issue_store <- new.env(parent = emptyenv())
    issue_store$values <- list()
    add_issue <- function(severity, code, message, n_rows = NA_integer_) {
        issue_store$values <- c(
            issue_store$values,
            list(.audit_issue(severity, code, message, n_rows))
        )
        invisible(NULL)
    }
    summary <- data.frame(
        metric = c("pairwise_rows", "stages_requested"),
        value = as.character(c(nrow(pairwise), length(stage_order))),
        stringsAsFactors = FALSE
    )
    coverage <- .empty_audit_coverage()
    replicate_cells <- data.frame(
        group = character(), stage = character(), n_samples = integer(),
        stringsAsFactors = FALSE
    )

    stage_order_valid <- length(stage_order) >= 3L &&
        !anyNA(stage_order) && !anyDuplicated(stage_order)
    if (!stage_order_valid) {
        add_issue(
            "error", "invalid_stage_order",
            "stage_order must contain at least three unique non-missing stages."
        )
    }
    stage_order <- as.character(stage_order)
    min_comparators_valid <- length(min_comparators) == 1L &&
        is.numeric(min_comparators) && !is.na(min_comparators) &&
        is.finite(min_comparators) && min_comparators >= 1L &&
        min_comparators == as.integer(min_comparators)
    if (!min_comparators_valid) {
        add_issue(
            "error", "invalid_min_comparators",
            "min_comparators must be one positive integer."
        )
    }
    min_comparators_value <- if (min_comparators_valid) {
        as.integer(min_comparators)
    } else {
        1L
    }

    required <- c(
        feature_col, gene_col, focal_col, comparator_col, stage_col,
        effect_col, q_col
    )
    optional <- c(gene_name_col, gene_q_col)
    missing_columns <- setdiff(
        c(required, optional[!is.null(optional)]), colnames(pairwise)
    )
    if (length(missing_columns)) {
        add_issue(
            "error", "missing_columns",
            paste(
                "Required columns are absent:",
                paste(missing_columns, collapse = ", ")
            ),
            length(missing_columns)
        )
        add_issue(
            "info", "manual_inference_check",
            paste(
                "Confirm focal-minus-comparator orientation and adjustment",
                "over the complete searched family; these cannot be audited."
            )
        )
        return(.finish_input_audit(
            summary, issue_store$values, coverage, replicate_cells
        ))
    }

    canonical <- data.frame(
        feature_id = as.character(pairwise[[feature_col]]),
        gene_id = as.character(pairwise[[gene_col]]),
        focal_group = as.character(pairwise[[focal_col]]),
        comparator_group = as.character(pairwise[[comparator_col]]),
        stage = as.character(pairwise[[stage_col]]),
        effect = .as_audit_numeric(pairwise[[effect_col]]),
        q_value = .as_audit_numeric(pairwise[[q_col]]),
        stringsAsFactors = FALSE
    )
    invalid_id <- !stats::complete.cases(canonical[seq_len(5L)]) |
        !nzchar(canonical$feature_id) | !nzchar(canonical$gene_id) |
        !nzchar(canonical$focal_group) |
        !nzchar(canonical$comparator_group) | !nzchar(canonical$stage)
    if (any(invalid_id)) {
        add_issue(
            "error", "missing_identifiers",
            "Identifier fields contain missing or empty values.",
            sum(invalid_id)
        )
    }
    self_contrast <- canonical$focal_group == canonical$comparator_group
    self_contrast[is.na(self_contrast)] <- FALSE
    if (any(self_contrast)) {
        add_issue(
            "error", "self_contrasts",
            "Focal and comparator groups are identical.",
            sum(self_contrast)
        )
    }
    unknown_stage <- !canonical$stage %in% stage_order
    unknown_stage[is.na(unknown_stage)] <- TRUE
    if (any(unknown_stage)) {
        add_issue(
            "error", "unknown_stages",
            "Pairwise stages are absent from stage_order.",
            sum(unknown_stage)
        )
    }
    invalid_effect <- !is.finite(canonical$effect)
    if (any(invalid_effect)) {
        add_issue(
            "warning", "missing_effects",
            "Non-finite effects make their feature-stage cells ineligible.",
            sum(invalid_effect)
        )
    }
    missing_q <- !is.finite(canonical$q_value)
    if (any(missing_q)) {
        add_issue(
            "warning", "missing_q_values",
            "Missing q-values make their feature-stage cells ineligible.",
            sum(missing_q)
        )
    }
    out_of_range_q <- is.finite(canonical$q_value) &
        (canonical$q_value < 0 | canonical$q_value > 1)
    if (any(out_of_range_q)) {
        add_issue(
            "error", "out_of_range_q_values",
            "Finite q-values must lie in [0, 1].",
            sum(out_of_range_q)
        )
    }
    if (!is.null(gene_q_col)) {
        gene_q <- .as_audit_numeric(pairwise[[gene_q_col]])
        invalid_gene_q <- is.finite(gene_q) & (gene_q < 0 | gene_q > 1)
        if (any(invalid_gene_q)) {
            add_issue(
                "error", "out_of_range_gene_q_values",
                "Finite gene-level q-values must lie in [0, 1].",
                sum(invalid_gene_q)
            )
        }
    }
    atomic_key <- .stable_key(
        canonical$feature_id, canonical$focal_group,
        canonical$comparator_group, canonical$stage
    )
    duplicate_rows <- duplicated(atomic_key) | duplicated(
        atomic_key, fromLast = TRUE
    )
    if (any(duplicate_rows)) {
        add_issue(
            "error", "duplicate_contrasts",
            "Duplicate feature/focal/comparator/stage rows must be resolved.",
            sum(duplicate_rows)
        )
    }
    valid_mapping <- !invalid_id
    mapping <- split(
        canonical$gene_id[valid_mapping],
        canonical$feature_id[valid_mapping]
    )
    conflicting <- vapply(
        mapping, function(value) length(unique(value)) != 1L, logical(1)
    )
    if (any(conflicting)) {
        add_issue(
            "error", "conflicting_gene_mapping",
            "Features map to more than one gene identifier.",
            sum(conflicting)
        )
    }

    observed_groups <- sort(unique(c(
        canonical$focal_group[!invalid_id],
        canonical$comparator_group[!invalid_id]
    )))
    groups <- if (is.null(expected_groups)) {
        observed_groups
    } else {
        unique(as.character(expected_groups))
    }
    valid_groups <- !is.na(groups) & nzchar(groups)
    if (any(!valid_groups)) {
        add_issue(
            "error", "invalid_expected_groups",
            "expected_groups contains missing or empty values."
        )
    }
    groups <- groups[valid_groups]
    if (any(!observed_groups %in% groups)) {
        add_issue(
            "error", "omitted_expected_groups",
            "expected_groups omits at least one observed group."
        )
    }
    if (min_comparators_valid &&
            length(groups) < min_comparators_value + 1L) {
        add_issue(
            "error", "insufficient_groups",
            "The group count is smaller than min_comparators + 1."
        )
    }
    observed_focal_groups <- sort(unique(
        canonical$focal_group[!invalid_id]
    ))
    focal_groups <- if (is.null(expected_focal_groups)) {
        observed_focal_groups
    } else {
        unique(as.character(expected_focal_groups))
    }
    valid_focal_groups <- !is.na(focal_groups) & nzchar(focal_groups)
    invalid_focal_members <- valid_focal_groups &
        !focal_groups %in% groups
    if (any(!valid_focal_groups) || any(invalid_focal_members)) {
        add_issue(
            "error", "invalid_expected_focal_groups",
            "expected_focal_groups must be non-missing members of groups."
        )
    }
    focal_groups <- focal_groups[
        valid_focal_groups & !invalid_focal_members
    ]
    if (any(!observed_focal_groups %in% focal_groups)) {
        add_issue(
            "error", "omitted_focal_groups",
            "expected_focal_groups omits at least one observed focal group."
        )
    }

    usable <- !invalid_id & !self_contrast & !unknown_stage
    usable_rows <- canonical[usable, , drop = FALSE]
    features <- sort(unique(usable_rows$feature_id))
    if (length(groups) && length(focal_groups) && stage_order_valid &&
            length(features) && min_comparators_valid) {
        coverage_grid <- expand.grid(
            focal_group = focal_groups,
            stage = stage_order,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        coverage_rows <- lapply(seq_len(nrow(coverage_grid)), function(index) {
            focal <- coverage_grid$focal_group[[index]]
            stage <- coverage_grid$stage[[index]]
            selected <- usable_rows[
                usable_rows$focal_group == focal &
                    usable_rows$stage == stage,
                , drop = FALSE
            ]
            by_feature <- split(selected$comparator_group, selected$feature_id)
            required_comparators <- setdiff(groups, focal)
            complete <- if (length(by_feature)) {
                vapply(by_feature, function(value) {
                    observed <- unique(value)
                    if (isTRUE(require_all_comparators)) {
                        setequal(observed, required_comparators)
                    } else {
                        length(observed) >= min_comparators_value
                    }
                }, logical(1))
            } else {
                logical()
            }
            expected_atomic <- length(features) * length(required_comparators)
            data.frame(
                focal_group = focal,
                stage = stage,
                expected_feature_cells = length(features),
                observed_feature_cells = length(by_feature),
                complete_feature_cells = sum(complete),
                completeness_fraction = sum(complete) / length(features),
                observed_atomic_rows = length(unique(.stable_key(
                    selected$feature_id, selected$comparator_group
                ))),
                expected_atomic_rows = expected_atomic,
                stringsAsFactors = FALSE
            )
        })
        coverage <- do.call(rbind, coverage_rows)
        incomplete_cells <- sum(
            coverage$expected_feature_cells - coverage$complete_feature_cells
        )
        if (incomplete_cells > 0L) {
            add_issue(
                "warning", "incomplete_contrast_family",
                paste(
                    "Some expected feature/focal/stage cells lack the",
                    "required directed comparator family."
                ),
                incomplete_cells
            )
        }
    }

    summary <- rbind(
        summary,
        data.frame(
            metric = c(
                "features", "genes", "groups_observed",
                "focal_groups_scanned"
            ),
            value = as.character(c(
                length(unique(canonical$feature_id[!invalid_id])),
                length(unique(canonical$gene_id[!invalid_id])),
                length(observed_groups),
                length(observed_focal_groups)
            )),
            stringsAsFactors = FALSE
        )
    )

    if (!is.null(usage)) {
        extracted <- tryCatch(
            .extract_usage(
                usage, col_data, assay_name, sample_col,
                group_col, sample_stage_col
            ),
            error = identity
        )
        if (inherits(extracted, "error")) {
            add_issue(
                "error", "invalid_usage_data", conditionMessage(extracted)
            )
        } else {
            metadata <- extracted$col_data
            matrix_value <- extracted$usage
            sample_groups <- as.character(metadata[[group_col]])
            sample_stages <- as.character(metadata[[sample_stage_col]])
            if (length(groups) && stage_order_valid) {
                cell_grid <- expand.grid(
                    group = groups, stage = stage_order,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
                )
                cell_grid$n_samples <- vapply(
                    seq_len(nrow(cell_grid)),
                    function(index) sum(
                        sample_groups == cell_grid$group[[index]] &
                            sample_stages == cell_grid$stage[[index]]
                    ),
                    integer(1)
                )
                replicate_cells <- cell_grid
                if (any(cell_grid$n_samples == 0L)) {
                    add_issue(
                        "error", "missing_replicate_cells",
                        "Usage metadata lacks expected group-stage cells.",
                        sum(cell_grid$n_samples == 0L)
                    )
                }
                if (any(cell_grid$n_samples == 1L)) {
                    add_issue(
                        "warning", "single_sample_cells",
                        paste(
                            "Some group-stage cells have one sample;",
                            "replicate separation is descriptive only."
                        ),
                        sum(cell_grid$n_samples == 1L)
                    )
                }
            }
            missing_features <- setdiff(features, rownames(matrix_value))
            if (length(missing_features)) {
                add_issue(
                    "warning", "features_missing_from_usage",
                    "Detected features absent from usage cannot be checked.",
                    length(missing_features)
                )
            }
            nonfinite_usage <- sum(!is.finite(matrix_value))
            if (nonfinite_usage) {
                add_issue(
                    "warning", "missing_usage_values",
                    paste(
                        "Non-finite usage values require an explicit",
                        "missing_values policy."
                    ),
                    nonfinite_usage
                )
            }
            summary <- rbind(
                summary,
                data.frame(
                    metric = c("usage_features", "usage_samples"),
                    value = as.character(dim(matrix_value)),
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    add_issue(
        "info", "manual_inference_check",
        paste(
            "Confirm focal-minus-comparator orientation and adjustment over",
            "the complete searched family; these cannot be audited."
        )
    )
    .finish_input_audit(
        summary, issue_store$values, coverage, replicate_cells
    )
}

#' Print a transient-DTU input audit
#'
#' @param x A `TransientDTUAudit` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.TransientDTUAudit <- function(x, ...) {
    cat("Transient DTU input audit\n")
    cat("status:", toupper(x$status), "\n")
    cat("issues:", nrow(x$issues), "\n")
    if (nrow(x$issues)) {
        counts <- table(as.character(x$issues$severity))
        cat(
            paste(names(counts), as.integer(counts), sep = "="),
            sep = "  "
        )
        cat("\n")
        print(as.data.frame(x$issues), row.names = FALSE)
    }
    invisible(x)
}
