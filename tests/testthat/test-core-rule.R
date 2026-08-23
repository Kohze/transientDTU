test_that("the bounded reciprocal fixture is detected exactly", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture)
    stage <- makeStageDTU(
        pairwise,
        fixture$stages,
        gene_name_col = "gene_name",
        gene_q_col = "gene_q"
    )
    episodes <- detectEpisodes(stage, gene_q_threshold = 0.05)

    expect_s4_class(stage, "DFrame")
    expect_s4_class(episodes, "DFrame")
    expect_equal(nrow(episodes), 2L)
    expect_setequal(episodes$direction, c("higher", "lower"))
    expect_true(all(episodes$start_stage == "s3"))
    expect_true(all(episodes$end_stage == "s3"))

    checked <- checkReplicateSeparation(
        episodes, fixture$usage, fixture$col_data, keep = "consistent"
    )
    expect_equal(nrow(checked), 2L)
    expect_true(all(checked$replicate_separation > 0))
    reciprocal <- annotateReciprocal(checked)
    expect_true(all(reciprocal$reciprocal_exchange))

    ranked <- rankCandidates(
        reciprocal,
        n = 6,
        replicate_consistent_only = TRUE
    )
    expect_equal(nrow(ranked), 1L)
    expect_equal(ranked$gene_id, "g1")
    expect_equal(ranked$rank, 1L)
})

test_that("consecutive states collapse and boundary states are excluded", {
    consecutive <- make_usage_fixture(event_stages = c("s3", "s4"))
    stage <- makeStageDTU(
        pairwise_from_fixture(consecutive), consecutive$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    episodes <- detectEpisodes(stage)
    expect_equal(nrow(episodes), 2L)
    expect_true(all(episodes$start_stage == "s3"))
    expect_true(all(episodes$end_stage == "s4"))
    expect_true(all(episodes$n_stages == 2L))

    boundary <- make_usage_fixture(event_stages = "s1")
    boundary_stage <- makeStageDTU(
        pairwise_from_fixture(boundary), boundary$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    expect_equal(nrow(detectEpisodes(boundary_stage)), 0L)
})

test_that("strict and inclusive thresholds follow the specification", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture, q_signal = 0.05)
    stage <- makeStageDTU(
        pairwise, fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    expect_equal(nrow(detectEpisodes(stage, q_threshold = 0.05)), 0L)

    pairwise$q_value[pairwise$q_value == 0.05] <- 0.001
    stage <- makeStageDTU(
        pairwise, fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    expect_equal(nrow(detectEpisodes(stage, effect_threshold = 0.30)), 2L)
    expect_equal(nrow(detectEpisodes(stage, effect_threshold = 0.31)), 0L)
})

test_that("missing comparison and flank evidence cannot create an episode", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture)
    pairwise <- pairwise[!(
        pairwise$feature_id == "g1_tx1" &
            pairwise$focal_group == "A" &
            pairwise$comparator_group == "C" &
            pairwise$stage == "s3"
    ), ]
    stage <- makeStageDTU(
        pairwise, fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    episodes <- detectEpisodes(stage)
    expect_false(any(episodes$feature_id == "g1_tx1"))

    pairwise <- pairwise_from_fixture(fixture)
    pairwise <- pairwise[!(
        pairwise$feature_id == "g1_tx1" &
            pairwise$focal_group == "A" & pairwise$stage == "s2"
    ), ]
    stage <- makeStageDTU(
        pairwise, fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    episodes <- detectEpisodes(stage)
    expect_false(any(episodes$feature_id == "g1_tx1"))
})

test_that("invalid and duplicate inputs fail early", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture)
    expect_error(
        makeStageDTU(rbind(pairwise, pairwise[1, ]), fixture$stages),
        "Duplicate"
    )
    pairwise$q_value[[1]] <- 1.2
    expect_error(makeStageDTU(pairwise, fixture$stages), "\\[0, 1\\]")

    conflicting_mapping <- pairwise_from_fixture(fixture)
    conflicting_mapping$gene_id[[1L]] <- "different_gene"
    expect_error(
        makeStageDTU(conflicting_mapping, fixture$stages),
        "exactly one gene"
    )
})
