test_that("complete workflow returns a valid interoperable result", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture)
    fit <- runTransientDTU(
        pairwise,
        fixture$stages,
        usage = fixture$usage,
        col_data = fixture$col_data,
        gene_name_col = "gene_name",
        gene_q_col = "gene_q",
        gene_q_threshold = 0.05
    )
    expect_s4_class(fit, "TransientDTUResult")
    expect_true(validObject(fit))
    expect_equal(nrow(episodeTable(fit)), 2L)
    expect_equal(nrow(candidateTable(fit)), 1L)
    expect_equal(decisionParameters(fit)$replicate_action, "filter")
    expect_true(all(c(
        "pairwise_rows", "episodes_retained", "candidate_groups"
    ) %in% diagnosticTable(fit)$metric))
})

test_that("replicate overlap filters episodes", {
    fixture <- make_usage_fixture()
    pairwise <- pairwise_from_fixture(fixture)
    fixture$usage["g1_tx1", fixture$col_data$group == "A" &
        fixture$col_data$stage == "s3"] <- c(0.39, 0.71)
    fixture$usage["g1_tx2", ] <- 1 - fixture$usage["g1_tx1", ]
    fit <- runTransientDTU(
        pairwise, fixture$stages,
        usage = fixture$usage, col_data = fixture$col_data,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    expect_equal(nrow(episodeTable(fit)), 0L)
})

test_that("incomplete replicate handling is explicit", {
    fixture <- make_usage_fixture(missing_sample = TRUE)
    pairwise <- pairwise_from_fixture(fixture)
    stage <- makeStageDTU(
        pairwise, fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    episodes <- detectEpisodes(stage)
    expect_error(
        checkReplicateSeparation(
            episodes, fixture$usage, fixture$col_data,
            incomplete = "error"
        )
    )
})

test_that("threshold grid reports stable selections", {
    fixture <- make_usage_fixture()
    stage <- makeStageDTU(
        pairwise_from_fixture(fixture), fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    )
    stability <- thresholdStability(
        stage,
        q_threshold = c(0.01, 0.05),
        effect_threshold = c(0.10, 0.20),
        comparator_tolerance = 0.10,
        flank_tolerance = 0.10,
        usage = fixture$usage,
        col_data = fixture$col_data,
        replicate_action = "filter"
    )
    expect_equal(nrow(stability$summary), 4L)
    expect_equal(stability$selection_frequency$gene_id, "g1")
    expect_equal(stability$selection_frequency$selection_frequency, 1)
})

test_that("SummarizedExperiment usage is accepted", {
    fixture <- make_usage_fixture()
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(usage = fixture$usage),
        colData = S4Vectors::DataFrame(fixture$col_data)
    )
    episodes <- detectEpisodes(makeStageDTU(
        pairwise_from_fixture(fixture), fixture$stages,
        gene_name_col = "gene_name", gene_q_col = "gene_q"
    ))
    checked <- checkReplicateSeparation(episodes, se, keep = "consistent")
    expect_equal(nrow(checked), 2L)
})
