context("import_data_path")

test_that("path is interpreted correctly", {
    f <- function(.data, path) {
        ## mock import_mtxzip() / import_loom()
        if (missing(.data) || missing(path)) {
            .import_data_path(.data, path, download = FALSE)
        } else {
            .import_data_path(.data, !!enquo(path), download = FALSE)
        }
    }

    expect_identical(f("foo"), "foo")
    expect_identical(f(path = "foo"), "foo")
    mypath <- "foo"
    expect_identical(f(mypath), "foo")

    df <- data.frame(foo = "bar")
    expect_identical(f(df, "foo"), "bar")
    expect_error(f(df), '"path" %in% names\\(.data\\) is not TRUE')

    ## tidy quotation
    expect_identical(f(df, foo), "bar")
    expect_identical(f(df, foo), "bar")

    ## default path look-up
    df <- data.frame(path = "bar")
    expect_identical(f(df), "bar")

    df <- data.frame(unknown = "bar")
    expect_error(f(df), '"path" %in% names\\(.data\\) is not TRUE')
})
