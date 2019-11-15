.coldata_n_elts <-
    function(x)
{
    vapply(x, function(column) length(unique(column)), integer(1))
}

#' Helper functions for simplifying colData
#'
#' @rdname colData_helpers
#'
#' @description `colData_tibble()` creates a tibble from `colData()`
#'
#' @param sce a `SummarizedExperiment` or derived object, e.g., a
#'     `SingleCellExperiment`.
#'
#' @return `colData_tibble()` returns a tibble representation of
#'     `colData()`.
#'
#' @export
colData_tibble <-
    function(sce)
{
    as_tibble(colData(sce))
}

#' @rdname colData_helpers
#'
#' @description `colData_column()` scans `colData()` columns to
#'     identify those with a single distinct value.
#'
#' @param abbreviate_colnames logical(1). When `TRUE`, abbreviate
#'     column names of `colData(sce)` using `names_abbreviate()`.
#'
#' @return `colData_column()` returns a two-column tibble. The first
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
colData_common <-
    function(sce, abbreviate_colnames = FALSE)
{
    stopifnot(
        is(sce, "SummarizedExperiment"),
        .is_scalar_logical(abbreviate_colnames)
    )

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    colData <- colData[n == 1L]
    value <- sapply(colData, unique)
    key <- names(colData)
    if (abbreviate_colnames)
        key <- names_abbreviate(setNames(nm = key))

    tibble(key, value)
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
    function(sce, abbreviate_colnames = FALSE)
{
    stopifnot(
        is(sce, "SummarizedExperiment"),
        .is_scalar_logical(colnames)
    )

    colData <- colData(sce)
    n <- .coldata_n_elts(colData)
    tbl <- as_tibble(colData[, n != 1L, drop = FALSE])
    if (abbreviate_colnames)
        names(tbl) <- names_abbreviate(tbl)

    tbl
}

#' @rdname colData_helpers
#'
#' @description `names_abbreviate()` abbreviates object `names()` to
#'     shortest common word suffix, using `sep` to define words.
#'
#' @param x Any object for which `names(x)` returns a character
#'     vector.
#'
#' @param sep character(1) regular expression identifying single
#'     characters to define words. The default splits a name into
#'     words at each occurrence of `_` or `.`.
#'
#' @param map logical(1) when FALSE (default) the return value is a
#'     vector of abbreviated names. When TRUE the return value is a
#'     tibble with columns `name` and `abbrev` corresponding to the
#'     original and abbreviated names.
#'
#' @return `names_abbreviate()` returns a vector of abbreviated names
#'     (`map = FALSE`, default) or a tibble (`map = TRUE`) with
#'     columns `name` and `abbrev` corresponding to the original and
#'     abbreviated names.
#'
#' @examples
#' nms <- setNames(nm = c("common_prefix.a_name", "common_prefix.a_title"))
#' names_abbreviate(nms)
#' names_abbreviate(nms, map=TRUE)
#'
#' @export
names_abbreviate <-
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
