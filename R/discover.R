.S3_BUCKET <-
    "https://s3.amazonaws.com/project-assets.data.humancellatlas.org/"

.s3_chr <-
    function(xml, elt)
{
    xpath <- paste0(
        "//d1:Key[text() != 'blacklist']/parent::d1:Contents/d1:",
        elt,
        "/text()"
    )
    nodeset <- xml_find_all(xml, xpath)
    as.character(nodeset)
}

#' Discover AWS S3 buckets with pre-computed HCA files
#'
#' @return A `tibble` describing available samples and the full path
#'     to the archive.
#'
#' @importFrom httr GET stop_for_status content
#' @importFrom xml2 xml_find_all
#' @importFrom tibble tibble
#'
#' @export
discover <-
    function()
{
    response <- GET(.S3_BUCKET)
    stop_for_status(response)

    xml <- content(response, encoding = "UTF-8")

    elt <- as.character(xml_find_all(xml, "//d1:IsTruncated/text()"))
    truncated <- !identical(elt, "false")
    ## FIXME: deal with truncated list

    key <- .s3_chr(xml, "Key")
    last_modified <- as.POSIXct(.s3_chr(xml, "LastModified"))
    etag <- .s3_chr(xml, "ETag")
    size <- as.numeric(.s3_chr(xml, "Size"))
    storage_class <- .s3_chr(xml, "StorageClass")

    tibble(
        file = basename(key),
        last_modified,
        etag = gsub('"', "", etag),
        size,
        storage_class,
        path = paste0(.S3_BUCKET, key)
    )
}
