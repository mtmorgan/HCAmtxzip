.data_n_elts <-
    function(x)
{
    vapply(x, function(column) length(unique(column)), integer(1))
}

#' @importFrom utils head tail
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
#' @rdname tbl_utilities
#'
#' @description `colTibble()` creates a tibble from `colData()`
#'
#' @param sce a `SummarizedExperiment` or derived object, e.g., a
#'     `SingleCellExperiment`.
#'
#' @return `colTibble()` returns a tibble representation of
#'     `colData()`.
#'
#' @export
colTibble <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    as_tibble(colData(sce))
}

#' @rdname tbl_utilities
#'
#' @return `rowTibble() returns a tibble representation of
#'     `rowData()`.
#'
#' @export
rowTibble <-
    function(sce)
{
    stopifnot(is(sce, "SummarizedExperiment"))

    as_tibble(rowData(sce))
}

#' @rdname tbl_utilities
#'
#' @description `constant()` scans `colData()` columns to
#'     identify those with a single distinct value.
#'
#' @param .data A `data.frame()` or `tibble()`.
#'
#' @return `constant()` returns a two-column tibble. The first column,
#'     `key` corresponds to `colData()` column names for which only a
#'     single distinct value is present. The second column `value` is
#'     the unique value.
#'
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom SingleCellExperiment rowData colData
#'
#' @export
constant <-
    function(.data)
{
    stopifnot(is.data.frame(.data))

    n <- .data_n_elts(.data)
    .data <- .data[n == 1L]
    value <- sapply(.data, unique)
    column <- names(.data)
    tibble(column, value)
}

#' @rdname tbl_utilities
#'
#' @description `brief()` produces a subset of `colData()`
#'     columns which contain more than one distinct value.
#'
#' @return `brief()` returns a tibble containing only those
#'     `colData()` columns with more than one value, and with column
#'     names abbreviated to the shortests common 'word' (using
#'     `[[:punct:]]+` as separators) suffixes.
#'
#' @export
brief <-
    function(.data)
{
    stopifnot(is.data.frame(.data))

    n <- .data_n_elts(.data)
    tbl <- as_tibble(.data[, n != 1L, drop = FALSE])
    names(tbl) <- .names_abbreviate(tbl)

    tbl
}
