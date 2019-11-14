.is_scalar_character <-
    function(x, allow.na = FALSE, allow.zchar = FALSE)
{
    is.character(x) && length(x) == 1L && (allow.na || !is.na(x)) &&
        (allow.zchar || nzchar(x))
}

.message <-
    function(...)
{
    ## FIXME: use futile.logger?
    message(...)
    TRUE
}

.MTX_ARCHIVE_FILES <- c(
    "features.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz", "cells.tsv.gz", 
    "barcodes.tsv.gz"
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

.download_zip <-
    function(path, destination)
{
    ## FIXME: use BiocFileCache
    if (!dir.exists(destination))
        dir.create(destination)

    destfile <- file.path(destination, basename(path))

    status <- download.file(path, destfile)
    if (!identical(status, 0L))
        stop("'download.file' failed with error code '", status, "'")

    destfile
}

.read_tsv <-
    function(path, ..., sep = "\t", header = FALSE)
{
    read.delim(path, ..., sep = sep, header = header)
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
#' @importFrom utils download.file read.delim unzip
#' @importFrom Matrix readMM
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
        path <- .download_zip(path, exdir)
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
    !verbose || .message("rowData / rowRanges")
    genes <- .read_tsv(ar[["genes"]], header=TRUE, row.names = "featurekey")
    genes <- makeGRangesFromDataFrame(
        genes,
        keep.extra.columns = TRUE,
        start.field = "featurestart", end.field = "featureend",
        starts.in.df.are.0based = TRUE
    )
    features <- .read_tsv(ar[["features"]]) # FIXME: add to `genes`?

    ## colData
    !verbose || .message("colData")
    cells <- .read_tsv(ar[["cells"]], header = TRUE, row.names = "cellkey")
    barcodes <-                         # FIXME: add to `cells`?
        .read_tsv(ar[["barcodes"]], blank.lines.skip = FALSE)


    ## assays
    !verbose || .message("assays")
    counts <- readMM(ar[["matrix"]])

    ## return value
    !verbose || .message("SingleCellExperiment")
    SingleCellExperiment(
        assays = list(counts = counts),
        colData = cbind(cells, barcodes),
        rowRanges = genes
    )
}


test <- c(
    "~/Downloads/Reprogrammed_Dendritic_Cells homo_sapiens 2019-11-08 16.12.mtx.zip",
    "https://s3.amazonaws.com/project-assets.data.humancellatlas.org/project-assets/project-matrices/116965f3-f094-4769-9d28-ae675c1b569c.homo_sapiens.mtx.zip"
)
