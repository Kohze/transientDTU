.assert_data_frame <- function(x, name) {
    if (!(is.data.frame(x) || methods::is(x, "DataFrame"))) {
        stop("'", name, "' must be a data.frame or DataFrame.", call. = FALSE)
    }
    as.data.frame(x, stringsAsFactors = FALSE)
}

.assert_columns <- function(x, columns, name) {
    missing <- setdiff(columns, colnames(x))
    if (length(missing)) {
        stop(
            "'", name, "' lacks required columns: ",
            paste(missing, collapse = ", "), ".",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.assert_number <- function(x, name, lower = -Inf, upper = Inf) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
            !is.finite(x) || x < lower || x > upper) {
        stop(
            "'", name, "' must be one finite number in [", lower, ", ",
            upper, "].",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.assert_count <- function(x, name, lower = 1L) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
            x != as.integer(x) || x < lower) {
        stop(
            "'", name, "' must be one integer >= ", lower, ".",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.assert_count_or_inf <- function(x, name, lower = 1L) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
            !(is.infinite(x) && x > 0) &&
            (!is.finite(x) || x != as.integer(x) || x < lower)) {
        stop(
            "'", name, "' must be one integer >= ", lower,
            " or positive Inf.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.as_id <- function(x, name, allow_na = FALSE) {
    value <- as.character(x)
    invalid <- is.na(value) | !nzchar(value)
    if (!allow_na && any(invalid)) {
        stop(
            "'", name, "' contains missing or empty identifiers.",
            call. = FALSE
        )
    }
    value
}

.stable_key <- function(...) {
    values <- list(...)
    do.call(paste, c(values, sep = "\034"))
}

.empty_episode_table <- function() {
    data.frame(
        feature_id = character(),
        gene_id = character(),
        gene_name = character(),
        focal_group = character(),
        start_stage = character(),
        end_stage = character(),
        start_index = integer(),
        end_index = integer(),
        n_stages = integer(),
        start_coordinate = numeric(),
        end_coordinate = numeric(),
        coordinate_span = numeric(),
        direction = character(),
        max_abs_usage_difference = numeric(),
        mean_abs_usage_difference = numeric(),
        weakest_component_q = numeric(),
        flanking_max_abs_difference = numeric(),
        n_required_flanks = integer(),
        n_finite_flanks = integer(),
        flanking_complete = logical(),
        gene_q = numeric(),
        replicate_separation = numeric(),
        replicate_consistent = logical(),
        reciprocal_exchange = logical(),
        stringsAsFactors = FALSE
    )
}

.empty_candidate_table <- function() {
    data.frame(
        rank = integer(),
        gene_id = character(),
        gene_name = character(),
        focal_group = character(),
        start_stage = character(),
        end_stage = character(),
        start_index = integer(),
        end_index = integer(),
        n_stages = integer(),
        start_coordinate = numeric(),
        end_coordinate = numeric(),
        coordinate_span = numeric(),
        n_features = integer(),
        feature_ids = character(),
        directions = character(),
        reciprocal_exchange = logical(),
        max_weakest_component_q = numeric(),
        max_abs_usage_difference = numeric(),
        min_replicate_separation = numeric(),
        stringsAsFactors = FALSE
    )
}

.order_na_last_desc <- function(x) {
    ifelse(is.na(x), -Inf, x)
}

.greater_equal_numeric <- function(x, threshold) {
    tolerance <- sqrt(.Machine$double.eps) *
        pmax(1, abs(x), abs(threshold))
    x > threshold | abs(x - threshold) <= tolerance
}
