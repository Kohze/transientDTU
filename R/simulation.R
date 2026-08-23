#' Simulate known transient and non-transient transcript-usage trajectories
#'
#' Generates two-part compositional usage for genes assigned to null,
#' monotonic, persistent-divergence, or bounded transient trajectories. Pairwise
#' evidence uses focal-minus-comparator mean differences and normal-reference
#' p-values under the known simulation noise, followed by a global adjustment.
#'
#' @param n_genes Number of two-feature genes.
#' @param n_stages Number of ordered stages; at least five.
#' @param groups Character group identifiers; at least three.
#' @param n_replicates Replicates per group-stage cell.
#' @param signal Absolute usage shift for the first feature in affected genes.
#' @param noise_sd Gaussian replicate noise on the usage scale.
#' @param transient_stage Interior stage index for bounded signals.
#' @param class_prob Named or ordered probabilities for `null`, `monotonic`,
#'   `persistent`, and `transient` genes.
#' @param missing_fraction Fraction of atomic pairwise rows removed completely
#'   at random after inference.
#' @param adjustment Multiplicity method passed to [stats::p.adjust()].
#' @param seed Integer random seed.
#'
#' @return A list containing `pairwise`, `usage`, `col_data`, `truth`,
#'   `gene_truth`, `stage_order`, and `groups`.
#'
#' @details
#' This simulation evaluates the post-inference decision layer. Its
#' normal-reference component p-values use the known generating noise and are
#' not intended to replace an RNA-seq count simulator or quantify transcript
#' assignment uncertainty.
#'
#' @examples
#' simulated <- simulateTransientDTU(n_genes = 20, seed = 11)
#' table(simulated$gene_truth$trajectory)
#' head(simulated$truth)
#'
#' @export
simulateTransientDTU <- function(
    n_genes = 100L,
    n_stages = 8L,
    groups = c("group1", "group2", "group3"),
    n_replicates = 2L,
    signal = 0.30,
    noise_sd = 0.015,
    transient_stage = ceiling(n_stages / 2) + 1L,
    class_prob = rep(0.25, 4L),
    missing_fraction = 0,
    adjustment = "BH",
    seed = 1L
) {
    .assert_count(n_genes, "n_genes", 4L)
    .assert_count(n_stages, "n_stages", 5L)
    .assert_count(n_replicates, "n_replicates", 2L)
    .assert_number(signal, "signal", 0.01, 0.45)
    .assert_number(noise_sd, "noise_sd", .Machine$double.eps, 0.25)
    .assert_count(transient_stage, "transient_stage", 2L)
    .assert_number(missing_fraction, "missing_fraction", 0, 0.9)
    groups <- unique(.as_id(groups, "groups"))
    if (length(groups) < 3L) {
        stop("At least three groups are required.", call. = FALSE)
    }
    if (transient_stage >= n_stages) {
        stop("'transient_stage' must have an earlier and a later flank.",
            call. = FALSE)
    }
    trajectories <- c("null", "monotonic", "persistent", "transient")
    if (is.null(names(class_prob))) names(class_prob) <- trajectories
    class_prob <- class_prob[trajectories]
    if (anyNA(class_prob) || any(class_prob < 0) || sum(class_prob) <= 0) {
        stop("'class_prob' must provide non-negative values for all classes.",
            call. = FALSE)
    }
    class_prob <- class_prob / sum(class_prob)
    if (!adjustment %in% stats::p.adjust.methods) {
        stop("Unknown multiplicity adjustment method.", call. = FALSE)
    }
    withr::with_seed(seed, {

    stage_order <- paste0("stage", seq_len(n_stages))
    col_data <- expand.grid(
        replicate = seq_len(n_replicates),
        group = groups,
        stage = stage_order,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    col_data$stage_index <- match(col_data$stage, stage_order)
    col_data$sample_id <- paste(
        col_data$group, col_data$stage, paste0("r", col_data$replicate),
        sep = "__"
    )
    rownames(col_data) <- col_data$sample_id

    gene_id <- sprintf("gene%04d", seq_len(n_genes))
    assignment <- sample(
        trajectories, n_genes, replace = TRUE, prob = class_prob
    )
    if (n_genes >= length(trajectories)) {
        assignment[seq_along(trajectories)] <- trajectories
    }
    focal_group <- rep(groups[[1L]], n_genes)
    base_usage <- stats::runif(n_genes, 0.35, 0.65)
    feature_id <- as.vector(rbind(
        paste0(gene_id, "_tx1"), paste0(gene_id, "_tx2")
    ))
    usage <- matrix(
        NA_real_, nrow = 2L * n_genes, ncol = nrow(col_data),
        dimnames = list(feature_id, col_data$sample_id)
    )

    for (gene_index in seq_len(n_genes)) {
        expected <- rep(base_usage[[gene_index]], nrow(col_data))
        is_focal <- col_data$group == focal_group[[gene_index]]
        stage_index <- col_data$stage_index
        trajectory <- assignment[[gene_index]]
        shift <- numeric(nrow(col_data))
        if (trajectory == "monotonic") {
            shift[is_focal] <- signal *
                (stage_index[is_focal] - 1) / (n_stages - 1)
        } else if (trajectory == "persistent") {
            shift[is_focal & stage_index >= transient_stage] <- signal
        } else if (trajectory == "transient") {
            shift[is_focal & stage_index == transient_stage] <- signal
        }
        first <- pmin(
            0.995,
            pmax(0.005, expected + shift + stats::rnorm(
                nrow(col_data), 0, noise_sd
            ))
        )
        first_id <- paste0(gene_id[[gene_index]], "_tx1")
        second_id <- paste0(gene_id[[gene_index]], "_tx2")
        usage[first_id, ] <- first
        usage[second_id, ] <- 1 - first
    }

    pair_rows <- vector(
        "list", 2L * n_genes * length(groups) *
            (length(groups) - 1L) * n_stages
    )
    row_index <- 0L
    known_se <- noise_sd * sqrt(2 / n_replicates)
    for (feature in rownames(usage)) {
        gene <- sub("_tx[12]$", "", feature)
        for (stage in stage_order) {
            for (focal in groups) {
                focal_values <- usage[
                    feature,
                    col_data$group == focal & col_data$stage == stage
                ]
                for (comparator in setdiff(groups, focal)) {
                    comparator_values <- usage[
                        feature,
                        col_data$group == comparator & col_data$stage == stage
                    ]
                    effect <- mean(focal_values) - mean(comparator_values)
                    p_value <- 2 * stats::pnorm(
                        -abs(effect / known_se)
                    )
                    row_index <- row_index + 1L
                    pair_rows[[row_index]] <- data.frame(
                        feature_id = feature,
                        gene_id = gene,
                        gene_name = gene,
                        focal_group = focal,
                        comparator_group = comparator,
                        stage = stage,
                        effect = effect,
                        p_value = p_value,
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
    pairwise <- do.call(rbind, pair_rows)
    pairwise$q_value <- stats::p.adjust(pairwise$p_value, method = adjustment)
    gene_min <- tapply(pairwise$p_value, pairwise$gene_id, min)
    gene_q <- stats::p.adjust(
        pmin(1, gene_min * n_stages * length(groups)),
        method = adjustment
    )
    pairwise$gene_q <- unname(gene_q[pairwise$gene_id])
    if (missing_fraction > 0) {
        keep <- stats::runif(nrow(pairwise)) >= missing_fraction
        pairwise <- pairwise[keep, , drop = FALSE]
    }
    rownames(pairwise) <- NULL

    transient_genes <- gene_id[assignment == "transient"]
    truth <- do.call(rbind, lapply(transient_genes, function(gene) {
        data.frame(
            feature_id = c(paste0(gene, "_tx1"), paste0(gene, "_tx2")),
            gene_id = gene,
            focal_group = groups[[1L]],
            start_stage = stage_order[[transient_stage]],
            end_stage = stage_order[[transient_stage]],
            direction = c("higher", "lower"),
            stringsAsFactors = FALSE
        )
    }))
    if (is.null(truth)) {
        truth <- data.frame(
            feature_id = character(), gene_id = character(),
            focal_group = character(), start_stage = character(),
            end_stage = character(), direction = character(),
            stringsAsFactors = FALSE
        )
    }
    gene_truth <- data.frame(
        gene_id = gene_id,
        trajectory = assignment,
        focal_group = focal_group,
        stringsAsFactors = FALSE
    )
    list(
        pairwise = pairwise,
        usage = usage,
        col_data = col_data,
        truth = truth,
        gene_truth = gene_truth,
        stage_order = stage_order,
        groups = groups,
        simulation_parameters = list(
            n_genes = n_genes,
            n_stages = n_stages,
            n_replicates = n_replicates,
            signal = signal,
            noise_sd = noise_sd,
            transient_stage = transient_stage,
            adjustment = adjustment,
            seed = seed
        )
    )
    })
}

.episode_identity <- function(x) {
    .stable_key(
        x$feature_id, x$gene_id, x$focal_group,
        x$start_stage, x$end_stage, x$direction
    )
}

#' Benchmark episode recovery under known truth
#'
#' @param n_simulations Number of independent simulations.
#' @param simulation_args Named arguments passed to [simulateTransientDTU()].
#' @param detector_args Named arguments passed to [runTransientDTU()].
#' @param seed First simulation seed; subsequent runs increment it by one.
#'
#' @return A list with per-run `results` and an aggregate `summary`, both as
#'   [S4Vectors::DataFrame] objects.
#'
#' @examples
#' benchmark <- benchmarkTransientDTU(
#'     n_simulations = 2,
#'     simulation_args = list(n_genes = 12),
#'     seed = 20
#' )
#' benchmark$summary
#'
#' @export
benchmarkTransientDTU <- function(
    n_simulations = 20L,
    simulation_args = list(),
    detector_args = list(),
    seed = 1L
) {
    .assert_count(n_simulations, "n_simulations")
    .assert_count(seed, "seed", 0L)
    if (!is.list(simulation_args) ||
        (length(simulation_args) > 0L &&
            (is.null(names(simulation_args)) ||
                any(!nzchar(names(simulation_args)))))) {
        stop("'simulation_args' must be a named list.", call. = FALSE)
    }
    if (!is.list(detector_args) ||
        (length(detector_args) > 0L &&
            (is.null(names(detector_args)) ||
                any(!nzchar(names(detector_args)))))) {
        stop("'detector_args' must be a named list.", call. = FALSE)
    }
    runs <- lapply(seq_len(n_simulations), function(index) {
        simulation_args$seed <- seed + index - 1L
        simulated <- do.call(simulateTransientDTU, simulation_args)
        fit_args <- utils::modifyList(
            list(
                pairwise = simulated$pairwise,
                stage_order = simulated$stage_order,
                usage = simulated$usage,
                col_data = simulated$col_data,
                panel_size = Inf
            ),
            detector_args
        )
        fit <- do.call(runTransientDTU, fit_args)
        detected <- as.data.frame(episodeTable(fit))
        truth_key <- .episode_identity(simulated$truth)
        detected_key <- .episode_identity(detected)
        true_positive <- sum(detected_key %in% truth_key)
        false_positive <- sum(!detected_key %in% truth_key)
        false_negative <- sum(!truth_key %in% detected_key)
        null_genes <- simulated$gene_truth$gene_id[
            simulated$gene_truth$trajectory == "null"
        ]
        false_null_genes <- unique(detected$gene_id[
            detected$gene_id %in% null_genes
        ])
        data.frame(
            simulation = index,
            seed = simulation_args$seed,
            truth_episodes = length(truth_key),
            detected_episodes = length(detected_key),
            true_positive = true_positive,
            false_positive = false_positive,
            false_negative = false_negative,
            precision = if (length(detected_key)) {
                true_positive / length(detected_key)
            } else {
                NA_real_
            },
            recall = if (length(truth_key)) {
                true_positive / length(truth_key)
            } else {
                NA_real_
            },
            any_false_positive = false_positive > 0L,
            null_gene_fpr = if (length(null_genes)) {
                length(false_null_genes) / length(null_genes)
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    })
    results <- do.call(rbind, runs)
    numeric_metrics <- c(
        "truth_episodes", "detected_episodes", "true_positive",
        "false_positive", "false_negative", "precision", "recall",
        "any_false_positive", "null_gene_fpr"
    )
    summary <- do.call(rbind, lapply(numeric_metrics, function(metric) {
        values <- as.numeric(results[[metric]])
        data.frame(
            metric = metric,
            mean = mean(values, na.rm = TRUE),
            median = stats::median(values, na.rm = TRUE),
            q025 = unname(stats::quantile(values, 0.025, na.rm = TRUE)),
            q975 = unname(stats::quantile(values, 0.975, na.rm = TRUE)),
            stringsAsFactors = FALSE
        )
    }))
    list(
        results = S4Vectors::DataFrame(results),
        summary = S4Vectors::DataFrame(summary)
    )
}
