context("colData_helpers")

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

    x <- setNames(nm = c("a_b_c", "d_b.c"))
    expect_equal(.names_abbreviate(x), c("a_b_c", "d_b.c"))

    x <- setNames(nm = c("a_b_c", "a_d.c"))
    expect_equal(.names_abbreviate(x), c("b_c", "d.c"))
})
