#' @importFrom BiocFileCache BiocFileCache bfcrpath
.download_file <-
    function(path)
{
    bfcrpath(BiocFileCache(), path)
}


.import_data_path <-
    function(.data, path, verbose = FALSE, download = TRUE)
{
    stopifnot(
        .is_scalar_logical(verbose),
        .is_scalar_logical(download)
    )

    if (!missing(.data) && missing(path)) {
        if (is.data.frame(.data)) {
            ## default for `path =`
            stopifnot("path" %in% names(.data))
            path <- pull(.data, "path")
        } else {
            path <- .data
        }
    } else if (missing(.data) && !missing(path)) {
        path <- path
    } else {
        path <- pull(.data, !!enquo(path))
    }

    path <- as.character(path)
    stopifnot(.is_scalar_character(path))

    ## download?
    if (download && startsWith(path, "http")) {
        !verbose || .message("download")
        path <- .download_file(path)
    }

    path
}
