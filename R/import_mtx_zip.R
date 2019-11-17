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

#' @importFrom BiocFileCache BiocFileCache bfcrpath
.download_zip <-
    function(path)
{
    bfcrpath(BiocFileCache(), path)
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
    v <- scan(
        path, list(integer(), integer(), numeric()), skip = 2,
        quiet = !verbose
    )
    sparseMatrix(v[[1]], v[[2]], x = v[[3]], dims = dims)
}

#' Import Human Cell Atlas '.mtx.zip' archives to SingleCellExperiment
#'
#' @rdname import_mtx_zip
#'
#' @param path character(1) path to remote (`http://` or `https://`)
#'     or local `.zip` archive.
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
#' url <- paste0(
#'     "https://s3.amazonaws.com/project-assets.data.humancellatlas.org/",
#'     "project-assets/project-matrices/",
#'     "116965f3-f094-4769-9d28-ae675c1b569c.homo_sapiens.mtx.zip"
#' )
#' sce <- import.mtxzip(url)
#' sce
#' }
#' @export
import.mtxzip <-
    function(path, ..., exdir = tempfile(), overwrite = FALSE, verbose = FALSE)
{
    stopifnot(
        .is_scalar_character(path), .is_scalar_character(exdir),
        startsWith(path, "http") || file.exists(path),
        !file.exists(exdir) || dir.exists(exdir),
        length(list(...)) == 0L
    )

    ## download?
    if (startsWith(path, "http")) {
        !verbose || .message("download")
        path <- .download_zip(path)
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


test <- c(
    "~/Downloads/Reprogrammed_Dendritic_Cells homo_sapiens 2019-11-08 16.12.mtx.zip",
    "https://s3.amazonaws.com/project-assets.data.humancellatlas.org/project-assets/project-matrices/116965f3-f094-4769-9d28-ae675c1b569c.homo_sapiens.mtx.zip"
)
