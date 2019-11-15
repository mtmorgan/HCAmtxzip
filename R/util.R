.is_scalar_character <-
    function(x, allow.na = FALSE, allow.zchar = FALSE)
{
    is.character(x) && length(x) == 1L && (allow.na || !is.na(x)) &&
        (allow.zchar || nzchar(x))
}

.is_scalar_integer <-
    function(x, allow.na = FALSE)
{
    is.integer(x) && length(x) == 1L && (allow.na || !is.na(x))
}

.is_scalar_logical <-
    function(x, allow.na = FALSE)
{
    is.logical(x) && length(x) == 1L && (allow.na || !is.na(x))
}

.message <-
    function(...)
{
    ## FIXME: use futile.logger?
    message(...)
    TRUE
}

#' @importFrom utils object.size

.object_size <-
    function(x)
{
    print(object.size(x), units = "auto")
}
