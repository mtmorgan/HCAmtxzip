.DSS_URL <- "https://dss.data.humancellatlas.org/v1"

.AZUL_SERVICE_URL <- "https://service.explore.data.humancellatlas.org"

.S3_BUCKET <-
    "https://s3.amazonaws.com/project-assets.data.humancellatlas.org/"

.AZUL_BUCKET <- "https://data.humancellatlas.org/"

.projects <-
    function(max = 1000L)
{
    stopifnot(.is_scalar_integer(max))
    query <- paste0(.AZUL_SERVICE_URL, "/repository/projects?size=", max)
    response <- GET(query)
    stop_for_status(response)

    pagination <- content(response)$pagination
    if (pagination$count != pagination$total)
        stop("[internal] too many projects available, contact maintainer")

    hits <- content(response)$hits
    entryId <- vapply(hits, `[[`, character(1), "entryId")
    projectTitle <-
        vapply(hits, function(x) x$projects[[1]]$projectTitle, character(1))
    tibble(projectTitle, entryId, hits)
}

#' @importFrom dplyr bind_cols bind_rows
#' @importFrom purrr discard keep modify_if
#' @importFrom jsonlite fromJSON
.file <-
    function(uuid, verbose = FALSE)
{
    !verbose || .message(uuid)
    query <- paste0(.DSS_URL, "/files/", uuid, "?replica=aws")
    path <- bfcrpath(BiocFileCache(), query)
    content <- fromJSON(path)

    kv = discard(content, is.data.frame) %>% unlist() %>% bind_rows()
    tbls = keep(kv, is.data.frame) %>% modify_if(is.data.frame, as_tibble)
    tbls <- Map(function(tbl, nm) {
        setNames(tibble(list(tbl)), nm)
    }, tbls, names(tbls))
    do.call(bind_cols, c(list(kv), tbls))
}

.files <-
    function(uuids, verbose = FALSE)
{
    responses <- lapply(uuids, .file, verbose = verbose)
    do.call(bind_rows, responses)
}

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

.buckets <-
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

    file <- basename(key)
    re <- "^([^.]+).(.*)$"

    tibble(
        entryId = sub(re, "\\1", file),
        fileFormat = sub(re, "\\2", file),
        size,
        path = paste0(.AZUL_BUCKET, key),
    )
}

#' Discover projects with pre-computed HCA files
#'
#' @return A `tibble` describing available projects and the full path
#'     to the archive.
#'
#' @examples
#' dd <- discover()
#' dd
#' dd %>% select(projectTitle)
#' dd %>%
#'     filter(grepl("^A single-cell reference", projectTitle)) %>%
#'     t()
#' path <- dd %>%
#'     filter(row_number() == which.min(size)) %>%
#'     pull(path)
#' sce <- import.mtxzip(path)
#'
#' @importFrom httr GET stop_for_status content
#' @importFrom xml2 xml_find_all
#' @importFrom dplyr "%>%" filter left_join
#' @importFrom tibble tibble as_tibble
#'
#' @export
discover <-
    function()
{
    projects <- .projects()
    buckets <- .buckets()
    suppressMessages(left_join(projects, buckets)) %>%
        filter(fileFormat %in% "mtx.zip") %>%
        select(-fileFormat)
}
