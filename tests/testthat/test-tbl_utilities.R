context("tbl_utilities")

test_that(".names_abbreviate works", {
    x <- setNames(nm = character())
    expect_equal(.names_abbreviate(x), character(0))

    x <- setNames(nm = NA_character_)
    expect_equal(.names_abbreviate(x), NA_character_)

    x <- setNames(nm = c("a", NA_character_))
    expect_equal(.names_abbreviate(x), c("a", NA_character_))

    x <- "a"
    expect_error(.names_abbreviate(x), "!is.null\\(names\\(x\\)\\) is not TRUE")

    x <- setNames(nm = c("b_a", "c_a"))
    expect_equal(.names_abbreviate(x), c("b_a", "c_a"))

    x <- setNames(nm = c("a_b", "a_c"))
    expect_equal(.names_abbreviate(x), c("b", "c"))

    x <- setNames(nm = c("a_b", "a.c"))
    expect_equal(.names_abbreviate(x), c("b", "c"))

    x <- setNames(nm = c("a_b_c", "d_b_c"))
    expect_equal(.names_abbreviate(x), c("a_b_c", "d_b_c"))

    ## '_' and '.' are word separators (from '[:punct:]')
    x <- setNames(nm = c("a_b_c", "a_d.c"))
    expect_equal(.names_abbreviate(x), c("b_c", "d.c"))

    ## '_' and '.' are treated differently in terms of matches
    x <- setNames(nm = c("a_b_c", "a_b.c", "a_b__c"))
    expect_equal(.names_abbreviate(x), c("b_c", "b.c", "b__c"))

    ## trailing punctuation
    x <- setNames(nm = c("b_c", "b_c_"))
    expect_equal(.names_abbreviate(x), c("c", "c_"))

    ## suffix of one matches all of another
    x <- setNames(nm = c("a.b.c", "b.c"))
    expect_equal(.names_abbreviate(x), unname(x))
})

test_that("constant() works", {
    df <- data.frame(x=1:2, y = 3)
    expected <- tibble(column = "y", value = 3)
    expect_identical(constant(df), expected)

    df <- data.frame(x=1:2, y = 3, z = "four", stringsAsFactors = FALSE)
    expected <- tibble(
        column = c("y", "z"),
        value = c(3, "four")
    )
    expect_identical(constant(df), expected)

    df <- data.frame(x=1:2)
    expected <- tibble(
        column = character(0),
        value = setNames(list(), character())
    )
    expect_identical(constant(df), expected)
})

test_that("brief() works", {
    tbl <- tibble(x=1:2, y = 3)
    expect_identical(brief(tbl), tbl["x"])

    tbl <- tibble(x=rep(2, 3), y = 3)
    expect_identical(brief(tbl), tbl[,FALSE])
})
