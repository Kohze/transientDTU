test_that("simulation is reproducible and compositional", {
    withr::local_seed(999)
    rng_before <- .Random.seed
    first <- simulateTransientDTU(n_genes = 12, seed = 100)
    expect_identical(.Random.seed, rng_before)
    second <- simulateTransientDTU(n_genes = 12, seed = 100)
    expect_identical(first$pairwise, second$pairwise)
    expect_equal(first$usage, second$usage)
    pairs <- matrix(first$usage, nrow = 2)
    expect_equal(colSums(pairs), rep(1, ncol(pairs)), tolerance = 1e-12)
    expect_setequal(
        unique(first$gene_truth$trajectory),
        c("null", "monotonic", "persistent", "transient")
    )
})

test_that("high-signal simulation recovers bounded truth", {
    simulated <- simulateTransientDTU(
        n_genes = 20, signal = 0.35, noise_sd = 0.005, seed = 101
    )
    fit <- runTransientDTU(
        simulated$pairwise,
        simulated$stage_order,
        usage = simulated$usage,
        col_data = simulated$col_data,
        gene_name_col = "gene_name",
        gene_q_col = "gene_q",
        panel_size = Inf
    )
    detected <- as.data.frame(episodeTable(fit))
    truth_keys <- transientDTU:::.episode_identity(simulated$truth)
    detected_keys <- transientDTU:::.episode_identity(detected)
    expect_true(all(truth_keys %in% detected_keys))
})

test_that("benchmark returns defined operating characteristics", {
    benchmark <- benchmarkTransientDTU(
        n_simulations = 2,
        simulation_args = list(
            n_genes = 12, signal = 0.35, noise_sd = 0.005
        ),
        detector_args = list(
            gene_name_col = "gene_name", gene_q_col = "gene_q"
        ),
        seed = 200
    )
    expect_s4_class(benchmark$results, "DFrame")
    expect_s4_class(benchmark$summary, "DFrame")
    expect_true(all(c("precision", "recall") %in% benchmark$summary$metric))

    empty_detector_args <- benchmarkTransientDTU(
        n_simulations = 1,
        simulation_args = list(n_genes = 8),
        detector_args = list(),
        seed = 201
    )
    expect_equal(nrow(empty_detector_args$results), 1L)
})

test_that("limma adapter returns the generic contract", {
    skip_if_not_installed("limma")
    simulated <- simulateTransientDTU(n_genes = 8, seed = 300)
    feature_data <- data.frame(
        gene_id = sub("_tx[12]$", "", rownames(simulated$usage)),
        row.names = rownames(simulated$usage)
    )
    pairwise <- fitPairwiseUsage(
        simulated$usage,
        simulated$col_data,
        feature_data,
        stage_order = simulated$stage_order,
        robust = FALSE
    )
    expect_true(all(c(
        "feature_id", "gene_id", "focal_group", "comparator_group",
        "stage", "effect", "q_value", "gene_q"
    ) %in% colnames(pairwise)))
    expect_equal(
        nrow(pairwise),
        nrow(simulated$usage) * length(simulated$stage_order) * 6L
    )
})
