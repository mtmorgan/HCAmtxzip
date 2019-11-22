.MTX_ARCHIVE_FILES <- c(
    "features.tsv.gz", "matrix.mtx.gz", "cells.tsv.gz", "barcodes.tsv.gz"
)

.is_mtx_archive <-
    function(x)
{
    msg <- character()

    test <- .MTX_ARCHIVE_FILES %in% basename(dir(x))
    if (!all(test)) {
        missing <- .MTX_ARCHIVE_FILES[!test]
        msg <- c(msg, "missing files:", paste0("  ", missing))
    }

    if (length(msg)) msg else TRUE
}

.read_tsv <-
    function(path, ..., stringsAsFactors = FALSE, sep = "\t", header = FALSE)
{
    read.delim(
        path, ...,
        stringsAsFactors = stringsAsFactors,
        sep = sep, header = header
    )
}

#' @importFrom Matrix sparseMatrix
.read_mtx <-
    function(path, verbose = FALSE)
{
    headers <- readLines(path, 2L)
    dims <- as.integer(strsplit(headers[2], " ")[[1]][1:2])
    !verbose || .message("dim: ", dims[1], " ", dims[2])
    v <- scan(
        path, list(integer(), integer(), numeric()), skip = 2,
        quiet = !verbose
    )
    sparseMatrix(v[[1]], v[[2]], x = v[[3]], dims = dims)
}

#' Import Human Cell Atlas '.mtx.zip' or '.loom' archives
#'
#' @rdname import
#'
#' @param .data (optional) When present, a data.frame or derived
#'     class, e.g., tibble, containing a single row, and a column
#'     specified by the `path` argument and containing a path to a
#'     remote (`http://` or `https://`) or local `.zip` archive.
#'
#'     `.data` can also be a `character(1)` argument used to provide
#'     the path to the `.zip` archive directly.
#'
#' @param path character(1) the path to the remote (`http://` or
#'     `https://`) or local `.zip` archive or, when `.data` is
#'     present, the column name (default `"path"`) in which the path
#'     is found.
#'
#' @param ... additional arguments, not supported.
#'
#' @param exdir character(1) directory in which to extract .zip archive.
#'
#' @param overwrite logical(1) overwrite existing files in `exdir`?
#'
#' @param verbose logical(1) report progress using `message()`.
#'
#' @return `SingleCellExperiment()` representing the expreession
#'     data. Rows represent features and columns cells; counts are
#'     represented in a sparse matrix.
#'
#' @importFrom utils read.delim unzip
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' mtx <-
#'     available("mtx.zip") %>%
#'     filter(size == min(size)) %>%
#'     import_mtxzip()
#' }
#' @export
import_mtxzip <-
    function(.data, path, ...,
             exdir = tempfile(), overwrite = FALSE, verbose = FALSE)
{
    stopifnot(
        length(list(...)) == 0L,
        .is_scalar_character(exdir),
        .is_scalar_logical(overwrite),
        overwrite || !file.exists(exdir),
        .is_scalar_logical(verbose)
    )

    if (missing(.data) || missing(path)) {
        path <- .import_data_path(.data, path, verbose = verbose)
    } else {
        path <- .import_data_path(.data, !!enquo(path), verbose = verbose)
    }

    ## unzip
    if (endsWith(path, ".zip")) {
        !verbose || .message("unzip")
        unzip(path, exdir = exdir, overwrite = overwrite, junkpaths = TRUE)
        path <- exdir
    }

    valid <- .is_mtx_archive(path)
    if (!isTRUE(valid))
        stop("invalid 'mtx' archive:\n  ", paste(valid, collapse = "\n  "))

    ar <- file.path(exdir, .MTX_ARCHIVE_FILES)
    names(ar) <- sub(".(tsv|mtx).gz", "", .MTX_ARCHIVE_FILES)

    ## rowData / rowRanges
    !verbose || .message("rowData")
    features <- .read_tsv(ar[["features"]], row.names = 1)

    ## colData
    !verbose || .message("colData")
    cells <- .read_tsv(ar[["cells"]], header = TRUE, row.names = "cellkey")
    if (!"barcode" %in% names(cells)) {
        barcodes <- .read_tsv(
            ar[["barcodes"]], blank.lines.skip = FALSE, col.names = "barcode"
        )
        cells <- cbind(cells, barcodes)
    }

    ## assays
    !verbose || .message("assays")
    counts <- .read_mtx(ar[["matrix"]], verbose)

    ## return value
    !verbose || .message("SingleCellExperiment")
    SingleCellExperiment(
        assays = list(counts = counts),
        colData = cells,
        rowData = features
    )
}
