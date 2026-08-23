make_usage_fixture <- function(
    event_stages = "s3",
    event_group = "A",
    missing_sample = FALSE,
    stages = paste0("s", 1:5)
) {
    groups <- c("A", "B", "C")
    col_data <- expand.grid(
        replicate = 1:2,
        group = groups,
        stage = stages,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    col_data$sample_id <- with(
        col_data, paste(group, stage, replicate, sep = "__")
    )
    rownames(col_data) <- col_data$sample_id
    usage <- matrix(
        NA_real_, nrow = 2, ncol = nrow(col_data),
        dimnames = list(c("g1_tx1", "g1_tx2"), col_data$sample_id)
    )
    first <- rep(0.4, nrow(col_data))
    first[col_data$group == event_group & col_data$stage %in% event_stages] <- 0.7
    first <- first + rep(c(-0.005, 0.005), length.out = length(first))
    usage["g1_tx1", ] <- first
    usage["g1_tx2", ] <- 1 - first
    if (missing_sample) {
        missing_column <- which(
            col_data$group == event_group &
                col_data$stage == event_stages[[1L]]
        )[[1L]]
        usage[1, missing_column] <- NA_real_
    }
    list(
        usage = usage,
        col_data = col_data,
        stages = stages,
        groups = groups
    )
}

pairwise_from_fixture <- function(fixture, q_signal = 0.001, q_null = 0.5) {
    rows <- list()
    index <- 0L
    for (feature in rownames(fixture$usage)) {
        for (stage in fixture$stages) {
            for (focal in fixture$groups) {
                focal_mean <- mean(fixture$usage[
                    feature,
                    fixture$col_data$group == focal &
                        fixture$col_data$stage == stage
                ], na.rm = TRUE)
                for (comparator in setdiff(fixture$groups, focal)) {
                    comparator_mean <- mean(fixture$usage[
                        feature,
                        fixture$col_data$group == comparator &
                            fixture$col_data$stage == stage
                    ], na.rm = TRUE)
                    effect <- focal_mean - comparator_mean
                    index <- index + 1L
                    rows[[index]] <- data.frame(
                        feature_id = feature,
                        gene_id = "g1",
                        gene_name = "GeneOne",
                        focal_group = focal,
                        comparator_group = comparator,
                        stage = stage,
                        effect = effect,
                        q_value = if (abs(effect) >= 0.2) q_signal else q_null,
                        gene_q = 0.001,
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
    do.call(rbind, rows)
}
