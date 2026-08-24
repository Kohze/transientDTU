test_that("methylation-expression simulation is reproducible and compositional", {
    withr::local_seed(999)
    rng_before <- .Random.seed
    first <- simulateMethylationExpression(seed = 42)
    expect_identical(.Random.seed, rng_before)
    second <- simulateMethylationExpression(seed = 42)

    expect_identical(first, second)
    expect_equal(length(first$stage_order), 8L)
    expect_equal(length(unique(first$cpg$position)), 24L)
    expect_true(all(first$observations$methylation > 0))
    expect_true(all(first$observations$methylation < 1))

    fraction_sum <- stats::aggregate(
        transcript_fraction ~ sample_id,
        data = first$observations,
        FUN = sum
    )
    expect_equal(fraction_sum$transcript_fraction, rep(1, nrow(fraction_sum)))
})

test_that("methylation-expression panel returns auditable summaries", {
    simulated <- simulateMethylationExpression(seed = 43)
    output <- tempfile(fileext = ".png")
    grDevices::png(output, width = 1200, height = 1500, res = 150)
    panel <- plotMethylationExpression(
        simulated$observations,
        cpg = simulated$cpg,
        stage_order = simulated$stage_order,
        main = "Test panel"
    )
    grDevices::dev.off()

    expect_true(file.exists(output))
    expect_gt(file.info(output)$size, 0)
    expect_s4_class(panel$stage_summary, "DFrame")
    expect_s4_class(panel$correlations, "DFrame")
    expect_s4_class(panel$trajectories, "DFrame")
    expect_s4_class(panel$cpg_summary, "DFrame")
    expect_equal(nrow(panel$correlations), 2L)
    expect_equal(nrow(panel$stage_summary), 16L)
    expect_true(all(is.finite(panel$correlations$correlation)))
})

test_that("methylation-expression panel validates its input contract", {
    simulated <- simulateMethylationExpression(seed = 44)
    expect_error(
        plotMethylationExpression(
            simulated$observations[, -which(
                names(simulated$observations) == "methylation"
            )]
        ),
        "lacks required columns"
    )
    expect_error(
        plotMethylationExpression(
            simulated$observations,
            stage_order = simulated$stage_order[-1L]
        ),
        "omits observed stages"
    )
    expect_error(
        simulateMethylationExpression(transcripts = "one"),
        "two unique identifiers"
    )
})
