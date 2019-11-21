context("import_mtx_zip")

test_that("path is interpreted correctly", {
    f <- function(.data, path) {
        ## mock import.mtxzip() 
        if (!missing(.data) && !missing(path)) {
            .import_mtxzip_path(.data, !!enquo(path))
        } else {
            .import_mtxzip_path(.data, path)
        }
    }

    expect_identical(f("foo"), "foo")
    expect_identical(f(path = "foo"), "foo")

    df <- data.frame(foo = "bar")
    expect_identical(f(df, "foo"), "bar")
    expect_error(f(df), 'argument \"path\" is missing, with no default')

    ## tidy quotation
    expect_identical(f(df, foo), "bar")
})
