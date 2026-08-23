#!/usr/bin/env Rscript

# Optional local integration test. The paper data are intentionally not bundled
# with transientDTU. Set TRANSIENTDTU_PAPER_DIR to the manuscript root before
# running this script from an installed or development package session.

paper_dir <- Sys.getenv("TRANSIENTDTU_PAPER_DIR", unset = NA_character_)
if (is.na(paper_dir) || !nzchar(paper_dir)) {
    stop(
        "Set TRANSIENTDTU_PAPER_DIR to the manuscript root containing ",
        "dtu_analysis/data.",
        call. = FALSE
    )
}
paper_dir <- normalizePath(paper_dir, winslash = "/", mustWork = TRUE)
data_dir <- file.path(paper_dir, "dtu_analysis", "data")
table_dir <- file.path(paper_dir, "dtu_analysis", "tables")
required <- c(
    file.path(data_dir, "transient_regional_stage_evaluations_all.csv"),
    file.path(data_dir, "transient_regional_filtered_isoform_fractions.csv"),
    file.path(data_dir, "transient_regional_isoform_episodes.csv"),
    file.path(table_dir, "transient_regional_top_candidates.csv")
)
if (any(!file.exists(required))) {
    stop("The paper directory lacks one or more regression inputs.",
        call. = FALSE)
}

stage_order <- c(
    "10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"
)
stage_source <- utils::read.csv(required[[1L]], check.names = FALSE)
stage_data <- data.frame(
    feature_id = stage_source$isoform_id,
    gene_id = stage_source$gene_id,
    gene_name = stage_source$gene_name,
    focal_group = stage_source$region,
    stage = as.character(stage_source$stage),
    stage_index = stage_source$stage_index,
    n_comparators = 2L,
    complete_comparisons = TRUE,
    direction_coherent =
        sign(stage_source$target_vs_other_1) ==
            sign(stage_source$target_vs_other_2) &
        sign(stage_source$target_vs_other_1) != 0,
    direction_sign = ifelse(
        sign(stage_source$target_vs_other_1) ==
            sign(stage_source$target_vs_other_2),
        sign(stage_source$target_vs_other_1), 0
    ),
    mean_effect = stage_source$target_difference,
    min_abs_effect = pmin(
        abs(stage_source$target_vs_other_1),
        abs(stage_source$target_vs_other_2)
    ),
    max_abs_effect = pmax(
        abs(stage_source$target_vs_other_1),
        abs(stage_source$target_vs_other_2)
    ),
    comparator_spread = abs(stage_source$other_region_difference),
    worst_q = stage_source$max_pair_q,
    gene_q = stage_source$gene_q,
    stringsAsFactors = FALSE
)
stage_data <- S4Vectors::DataFrame(stage_data)
attr(stage_data, "stage_order") <- stage_order

episodes <- detectEpisodes(
    stage_data,
    q_threshold = 0.05,
    effect_threshold = 0.10,
    comparator_tolerance = 0.10,
    flank_tolerance = 0.10,
    gene_q_threshold = 0.05,
    min_comparators = 2L
)
cat("Episodes before replicate filtering: ", nrow(episodes), "\n", sep = "")
if (!all(c("NR_131127", "NR_131128") %in% episodes$feature_id)) {
    cat("C230072F16Rik episodes were absent before replicate filtering.\n")
}

fractions <- utils::read.csv(required[[2L]], check.names = FALSE)
usage <- as.matrix(fractions[, -1L, drop = FALSE])
storage.mode(usage) <- "numeric"
rownames(usage) <- fractions$isoform_id
parts <- strsplit(colnames(usage), "__", fixed = TRUE)
sample_data <- data.frame(
    sample_id = colnames(usage),
    group = vapply(parts, `[[`, character(1), 1L),
    stage_replicate = vapply(parts, `[[`, character(1), 2L),
    stringsAsFactors = FALSE
)
sample_data$stage <- sub("_[12]$", "", sample_data$stage_replicate)
rownames(sample_data) <- sample_data$sample_id

episodes <- checkReplicateSeparation(
    episodes,
    usage,
    sample_data,
    sample_col = "sample_id",
    group_col = "group",
    stage_col = "stage",
    incomplete = "drop",
    keep = "consistent"
)
episodes <- annotateReciprocal(episodes)
observed <- as.data.frame(episodes)
expected <- utils::read.csv(required[[3L]], check.names = FALSE)

identity_columns <- c(
    "feature_id", "gene_id", "gene_name", "focal_group", "start_stage",
    "end_stage", "n_stages", "direction", "replicate_consistent",
    "reciprocal_exchange"
)
expected_identity <- data.frame(
    feature_id = expected$isoform_id,
    gene_id = expected$gene_id,
    gene_name = expected$gene_name,
    focal_group = expected$region,
    start_stage = as.character(expected$start_stage),
    end_stage = as.character(expected$end_stage),
    n_stages = expected$n_stages,
    direction = expected$direction,
    replicate_consistent = expected$replicate_consistent,
    reciprocal_exchange = expected$reciprocal_exchange,
    stringsAsFactors = FALSE
)
sort_key <- function(value) {
    do.call(order, value[c(
        "gene_id", "feature_id", "focal_group", "start_stage",
        "end_stage", "direction"
    )])
}
observed <- observed[sort_key(observed), , drop = FALSE]
expected_identity <- expected_identity[
    sort_key(expected_identity), , drop = FALSE
]
expected <- expected[order(
    expected$gene_id, expected$isoform_id, expected$region,
    expected$start_stage, expected$end_stage, expected$direction
), , drop = FALSE]
rownames(observed) <- NULL
rownames(expected_identity) <- NULL
rownames(expected) <- NULL

cat(
    "Regression counts: detected=", nrow(observed),
    ", archived=", nrow(expected), "\n", sep = ""
)
observed_key <- do.call(paste, c(observed[identity_columns[1:6]], sep = "__"))
expected_key <- do.call(
    paste, c(expected_identity[identity_columns[1:6]], sep = "__")
)
if (!setequal(observed_key, expected_key)) {
    cat("Missing archived keys:\n")
    print(expected_key[!expected_key %in% observed_key])
    cat("Unexpected detected keys:\n")
    print(observed_key[!observed_key %in% expected_key])
}

stopifnot(
    nrow(observed) == 1348L,
    identical(observed[identity_columns], expected_identity[identity_columns]),
    isTRUE(all.equal(
        observed$max_abs_usage_difference,
        expected$max_abs_usage_difference,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        observed$mean_abs_usage_difference,
        expected$mean_abs_usage_difference,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        observed$weakest_component_q,
        expected$worst_pair_q,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        observed$flanking_max_abs_difference,
        expected$flanking_max_abs_difference,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        observed$replicate_separation,
        expected$replicate_separation,
        tolerance = 1e-12
    ))
)

ranked <- as.data.frame(rankCandidates(
    episodes,
    n = 6L,
    reciprocal_only = TRUE,
    replicate_consistent_only = TRUE,
    unique_genes = TRUE
))
expected_rank <- utils::read.csv(required[[4L]], check.names = FALSE)
stopifnot(
    identical(ranked$gene_id, expected_rank$gene_id),
    identical(ranked$gene_name, expected_rank$gene_name),
    identical(ranked$focal_group, expected_rank$region),
    identical(ranked$start_stage, as.character(expected_rank$start_stage)),
    identical(ranked$end_stage, as.character(expected_rank$end_stage)),
    identical(ranked$n_features, expected_rank$isoforms),
    isTRUE(all.equal(
        ranked$max_weakest_component_q,
        expected_rank$max_worst_pair_q,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        ranked$max_abs_usage_difference,
        expected_rank$max_abs_usage_difference,
        tolerance = 1e-12
    )),
    isTRUE(all.equal(
        ranked$min_replicate_separation,
        expected_rank$min_replicate_separation,
        tolerance = 1e-12
    ))
)

cat(
    "PASS: transientDTU reproduced 1,348 archived episodes and the exact ",
    "six-gene ranking.\n",
    sep = ""
)
