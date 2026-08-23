test_that("input audit separates ready, review, and fail states", {
    simulated <- simulateTransientDTU(n_genes = 8, seed = 401)
    ready <- auditTransientInput(
        simulated$pairwise,
        simulated$stage_order,
        usage = simulated$usage,
        col_data = simulated$col_data
    )
    expect_s3_class(ready, "TransientDTUAudit")
    expect_identical(ready$status, "ready")
    expect_true(ready$ready)
    expect_true(all(ready$coverage$completeness_fraction == 1))

    incomplete <- simulated$pairwise[-1L, , drop = FALSE]
    review <- auditTransientInput(incomplete, simulated$stage_order)
    expect_identical(review$status, "review")
    expect_true("incomplete_contrast_family" %in% review$issues$code)

    duplicated <- rbind(simulated$pairwise, simulated$pairwise[1L, ])
    failed <- auditTransientInput(duplicated, simulated$stage_order)
    expect_identical(failed$status, "fail")
    expect_true("duplicate_contrasts" %in% failed$issues$code)

    targeted <- simulated$pairwise[
        simulated$pairwise$focal_group == simulated$groups[[1L]],
        , drop = FALSE
    ]
    targeted_audit <- auditTransientInput(
        targeted,
        simulated$stage_order,
        expected_groups = simulated$groups,
        expected_focal_groups = simulated$groups[[1L]]
    )
    expect_identical(targeted_audit$status, "ready")
    expect_setequal(
        unique(targeted_audit$coverage$focal_group),
        simulated$groups[[1L]]
    )

    out_of_range <- simulated$pairwise
    out_of_range$q_value[[1L]] <- 2
    invalid_q <- auditTransientInput(out_of_range, simulated$stage_order)
    expect_identical(invalid_q$status, "fail")
    expect_true("out_of_range_q_values" %in% invalid_q$issues$code)

    invalid_control <- auditTransientInput(
        simulated$pairwise,
        c("duplicated", "duplicated"),
        min_comparators = "two"
    )
    expect_identical(invalid_control$status, "fail")
    expect_true(all(c(
        "invalid_stage_order", "invalid_min_comparators"
    ) %in% invalid_control$issues$code))
})

test_that("two-group designs are first-class with one comparator", {
    simulated <- simulateTransientDTU(
        n_genes = 8, groups = c("treated", "control"), seed = 402
    )
    audit <- auditTransientInput(
        simulated$pairwise,
        simulated$stage_order,
        min_comparators = 1L
    )
    expect_identical(audit$status, "ready")
    fit <- runTransientDTU(
        simulated$pairwise,
        simulated$stage_order,
        min_comparators = 1L,
        reciprocal_only = FALSE,
        panel_size = Inf
    )
    expect_s4_class(fit, "TransientDTUResult")
    expect_gt(nrow(episodeTable(fit)), 0L)
})

test_that("episode duration, flank width, and coordinates are explicit", {
    long_fixture <- make_usage_fixture(
        event_stages = "s4", stages = paste0("s", 1:7)
    )
    stage <- makeStageDTU(
        pairwise_from_fixture(long_fixture), long_fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    coordinates <- c(s1 = 0, s2 = 1, s3 = 2, s4 = 4, s5 = 8, s6 = 16, s7 = 32)
    episodes <- detectEpisodes(
        stage, flank_width = 2L, stage_coordinates = coordinates
    )
    expect_equal(nrow(episodes), 2L)
    expect_true(all(episodes$n_required_flanks == 4L))
    expect_true(all(episodes$start_coordinate == 4))
    expect_true(all(episodes$coordinate_span == 0))
    expect_equal(
        nrow(detectEpisodes(stage, min_episode_stages = 2L)),
        0L
    )
    expect_error(
        detectEpisodes(stage, stage_coordinates = rev(seq_along(coordinates))),
        "strictly increasing"
    )

    two_stage_fixture <- make_usage_fixture(event_stages = c("s3", "s4"))
    two_stage <- makeStageDTU(
        pairwise_from_fixture(two_stage_fixture), two_stage_fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    expect_equal(
        nrow(detectEpisodes(two_stage, max_episode_stages = 1L)),
        0L
    )
})

test_that("replicate summaries and missing-value policies are selectable", {
    fixture <- make_usage_fixture()
    episodes <- detectEpisodes(makeStageDTU(
        pairwise_from_fixture(fixture), fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    ))
    focal_columns <- fixture$col_data$group == "A" &
        fixture$col_data$stage == "s3"
    fixture$usage["g1_tx1", focal_columns] <- c(0.35, 0.70)
    fixture$usage["g1_tx2", ] <- 1 - fixture$usage["g1_tx1", ]
    complete <- checkReplicateSeparation(
        episodes, fixture$usage, fixture$col_data
    )
    median_check <- checkReplicateSeparation(
        episodes, fixture$usage, fixture$col_data,
        replicate_method = "median"
    )
    expect_false(any(complete$replicate_consistent))
    expect_true(all(median_check$replicate_consistent))

    missing <- make_usage_fixture(missing_sample = TRUE)
    missing_episodes <- detectEpisodes(makeStageDTU(
        pairwise_from_fixture(missing), missing$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    ))
    default_missing <- checkReplicateSeparation(
        missing_episodes, missing$usage, missing$col_data
    )
    omit_missing <- checkReplicateSeparation(
        missing_episodes, missing$usage, missing$col_data,
        missing_values = "omit", min_replicates = 1L
    )
    strict_count <- checkReplicateSeparation(
        missing_episodes, missing$usage, missing$col_data,
        missing_values = "omit", min_replicates = 2L
    )
    expect_true(anyNA(default_missing$replicate_separation))
    expect_false(anyNA(omit_missing$replicate_separation))
    expect_true(anyNA(strict_count$replicate_separation))
})

test_that("threshold stability accepts fixed episode controls", {
    fixture <- make_usage_fixture()
    stage <- makeStageDTU(
        pairwise_from_fixture(fixture), fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    stability <- thresholdStability(
        stage,
        q_threshold = 0.05,
        effect_threshold = 0.10,
        comparator_tolerance = 0.10,
        flank_tolerance = 0.10,
        episode_args = list(min_episode_stages = 2L)
    )
    expect_true(all(stability$summary$episodes == 0L))
    expect_error(
        thresholdStability(stage, episode_args = list(q_threshold = 0.01)),
        "Invalid 'episode_args'"
    )
})
