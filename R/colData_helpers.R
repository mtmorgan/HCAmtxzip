.coldata_n_elts <-
    function(x)
{
    vapply(x, function(column) length(unique(column)), integer(1))
}

#' Helper functions for simpilfying colData
#'
#' @rdname colData_helpers
#'
#' @description `colData_column()` scans `colData()` columns to
#'     identify those with a single distinct value.
#'
#' @return `colData_column()` returns a two-column tibble. The first
#'     column, `key` corresponds to `colData()` column names for which
#'     only a single distinct value is present. The second column
#'     `value` is the unique value.
#'
#' @export
colData_common <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    colData <- colData[n == 1L]
    value <- sapply(colData, unique)
    tibble(key = names(colData), value)
}

#' @rdname colData_helpers
#'
#' @description `colData_distinct()` produces a subset of `colData()`
#'     columns which contain more than one distinct value.
#'
#' @return `colData_distinct()` returns a tibble containing only those
#'     `colData()` columns with more than one value.
#'
#' @export
colData_distinct <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    as_tibble(colData[, n != 1L, drop = FALSE])
}        
