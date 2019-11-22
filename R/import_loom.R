#' @rdname import
#'
#' @examples
#' \dontrun{
#' loom <-
#'     available("loom") %>%
#'     filter(size == min(size)) %>%
#'     import_loom()
#' }
#' @export
import_loom <-
    function(.data, path, ..., verbose = FALSE)
{
    stopifnot(
        length(list(...)) == 0L,
        .is_scalar_logical(verbose)
    )

    if (missing(.data) || missing(path)) {
        path <- .import_data_path(.data, path, verbose = verbose)
    } else {
        path <- .import_data_path(.data, !!enquo(path), verbose = verbose)
    }

    !verbose || .message("LoomExperiment")
    LoomExperiment::import(LoomExperiment::LoomFile(path))
}
