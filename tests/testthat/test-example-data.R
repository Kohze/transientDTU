test_that("bundled example follows Bioconductor container conventions", {
    data("transientExampleSE", package = "transientDTU")
    data("transientExamplePairwise", package = "transientDTU")

    expect_s4_class(transientExampleSE, "SummarizedExperiment")
    expect_true(validObject(transientExampleSE))
    expect_identical(
        SummarizedExperiment::assayNames(transientExampleSE), "usage"
    )
    expect_equal(dim(transientExampleSE), c(24L, 48L))
    expect_identical(
        rownames(SummarizedExperiment::rowData(transientExampleSE)),
        rownames(transientExampleSE)
    )
    expect_identical(
        rownames(SummarizedExperiment::colData(transientExampleSE)),
        colnames(transientExampleSE)
    )
    expect_true(all(c(
        "feature_id", "gene_id", "gene_name"
    ) %in% names(SummarizedExperiment::rowData(transientExampleSE))))
    expect_true(all(c(
        "sample_id", "group", "stage", "replicate", "stage_index"
    ) %in% names(SummarizedExperiment::colData(transientExampleSE))))
    expect_s4_class(transientExamplePairwise, "DataFrame")
    expect_equal(nrow(transientExamplePairwise), 1152L)
    expect_false(anyNA(transientExamplePairwise$gene_q))
})

test_that("bundled example completes the documented workflow", {
    data("transientExampleSE", package = "transientDTU")
    data("transientExamplePairwise", package = "transientDTU")
    example_metadata <- S4Vectors::metadata(transientExampleSE)

    expect_true(all(c(
        "stage_order", "groups", "truth", "gene_truth",
        "simulation_parameters", "provenance"
    ) %in% names(example_metadata)))

    audit <- auditTransientInput(
        transientExamplePairwise,
        example_metadata$stage_order,
        usage = transientExampleSE
    )
    expect_identical(audit$status, "ready")

    fit <- runTransientDTU(
        transientExamplePairwise,
        stage_order = example_metadata$stage_order,
        usage = transientExampleSE,
        gene_name_col = "gene_name",
        gene_q_col = "gene_q",
        gene_q_threshold = 0.05,
        panel_size = 6
    )
    expect_s4_class(fit, "TransientDTUResult")
    expect_true(validObject(fit))
    expect_equal(nrow(episodeTable(fit)), 4L)
    expect_equal(candidateTable(fit)$gene_id, c("gene0004", "gene0005"))
})

test_that("legacy and Bioconductor examples describe the same simulation", {
    data("transientExample", package = "transientDTU")
    data("transientExampleSE", package = "transientDTU")
    data("transientExamplePairwise", package = "transientDTU")

    expect_equal(
        SummarizedExperiment::assay(transientExampleSE, "usage"),
        transientExample$usage
    )
    expect_equal(
        as.data.frame(transientExamplePairwise),
        transientExample$pairwise
    )
    expect_identical(
        S4Vectors::metadata(transientExampleSE)$stage_order,
        transientExample$stage_order
    )
})
