#' Useful functions for interacting with Human Cell Atlas resources
#'
#' @rdname hca_utilities
#'
#' @description `replica()` sets or returns the replica for data retrieval.
#'
#' @param x character(1) replica to use, either `"aws"` or `"gcp"`.
#'
#' @return character(1) the (updated) replica
#'
#' @examples
#' replica()
#'
#' @export
replica <- local({
    replica <- "aws"
    function(x = c("aws", "gcp")) {
        if (!missing(x))
            replica <<- match.arg(x)
        replica
    }
})
