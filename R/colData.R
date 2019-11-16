.coldata_n_elts <-
    function(x)
{
    vapply(x, function(column) length(unique(column)), integer(1))
}

.names_abbreviate <-
    function(x, sep = "[_\\.]", map = FALSE)
{
    stopifnot(
        !is.null(names(x)),
        .is_scalar_character(sep),
        .is_scalar_logical(map)
    )

    nms <- names(x)
    elts <- strsplit(nms, sep)

    repeat {
        lens <- lengths(elts)
        if (all(lens == 1L))
            break
        abbrev <- mapply(`[[`, elts, lens)
        tbl <- table(abbrev)
        unique <- abbrev %in% names(tbl)[tbl == 1L]
        elts[unique] <- as.list(abbrev[unique])
        elts[!unique] <- Map(function(x, len) {
            c(head(x, -2), paste(x[[len - 1L]], x[[len]], sep="."))
        }, elts[!unique], lens[!unique])
    }

    elts <- substring(nms, nchar(nms) - nchar(unlist(elts)) + 1L)
    if (map) {
        tibble(name = nms, abbrev = elts)
    } else elts
}

#' Helper functions for simplifying colData
#'
#' @rdname colData
#'
#' @description `colDataTibble()` creates a tibble from `colData()`
#'
#' @param sce a `SummarizedExperiment` or derived object, e.g., a
#'     `SingleCellExperiment`.
#'
#' @return `colDataTibble()` returns a tibble representation of
#'     `colData()`.
#'
#' @export
colDataTibble <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    as_tibble(colData(sce))
}

#' @rdname colData
#'
#' @description `colDataConstants()` scans `colData()` columns to
#'     identify those with a single distinct value.
#'
#' @return `colDataConstants()` returns a two-column tibble. The first
#'     column, `key` corresponds to `colData()` column names for which
#'     only a single distinct value is present. The second column
#'     `value` is the unique value.
#'
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom SingleCellExperiment colData
#'
#' @export
colDataConstants <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    colData <- colData[n == 1L]
    value <- sapply(colData, unique)
    column <- names(colData)
    tibble(column, value)
}

#' @rdname colData
#'
#' @description `colDataBrief()` produces a subset of `colData()`
#'     columns which contain more than one distinct value.
#'
#' @return `colDataBrief()` returns a tibble containing only those
#'     `colData()` columns with more than one value, and with column
#'     names abbreviated to the shortests common 'word' (using `_` or
#'     `.` as separators) suffixes.
#'
#' @export
colDataBrief <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    tbl <- as_tibble(colData[, n != 1L, drop = FALSE])
    names(tbl) <- .names_abbreviate(tbl)

    tbl
}
