.coldata_n_elts <-
    function(x)
{
    vapply(x, function(column) length(unique(column)), integer(1))
}

.names_abbreviate <-
    function(x, map = FALSE)
{
    stopifnot(
        !is.null(names(x)),
        .is_scalar_logical(map)
    )

    nms <- names(x)
    ## allow for trailing punctation -- only split if the next
    ## charcter is not from punct other than at end of line
    elts <- strsplit(nms, "[[:punct:]]+(?=[^[:punct:]])", perl = TRUE)
    seps <- regmatches(nms, gregexpr("[[:punct:]]+(?=[^$])", nms, perl=TRUE))

    repeat {
        lens <- lengths(elts)
        if (all(lens == 1L))
            break
        abbrev <- lapply(elts, tail, 1)
        tbl <- table(unlist(abbrev))
        unique <- unlist(abbrev) %in% names(tbl)[tbl == 1L]
        elts[unique] <- abbrev[unique]
        elts[!unique] <- Map(function(elt, sep) {
            c(head(elt, -2), paste0(tail(elt, 2), collapse=sep))
        }, elts[!unique], lapply(seps[!unique], tail, 1))
        seps <- lapply(seps, head, -1L)
    }

    elts <- as.character(unlist(elts)) # as.character for length(x) == 0
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
#' @return `rowDataTibble() returns a tibble representation of
#'     `rowData()`.
#'
#' @export
rowDataTibble <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    as_tibble(rowData(scde))
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
#'     names abbreviated to the shortests common 'word' (using
#'     `[[:punct:]]+` as separators) suffixes.
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
