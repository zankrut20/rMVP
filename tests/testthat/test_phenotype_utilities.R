# tests/testthat/test_phenotype_utilities.R
#
# Unit tests for Phase 1 phenotype utility functions:
#   detect_family(), recode_phenotype(), validate_binary_phenotype()
#
# Run with: devtools::test(filter = "phenotype_utilities")

library(testthat)
library(rMVP)

# ---------------------------------------------------------------------------
# detect_family() tests
# ---------------------------------------------------------------------------

test_that("detect_family identifies {0,1} as binomial", {
    expect_equal(detect_family(c(0, 1, 1, 0, 0), verbose = FALSE), "binomial")
})

test_that("detect_family identifies {1,2} as binomial", {
    expect_equal(detect_family(c(1, 2, 1, 2, 2), verbose = FALSE), "binomial")
})

test_that("detect_family identifies {-1,1} as binomial", {
    expect_equal(detect_family(c(-1, 1, -1, 1), verbose = FALSE), "binomial")
})

test_that("detect_family returns gaussian for continuous phenotype", {
    expect_equal(detect_family(c(1.2, 3.5, 2.8, 4.1), verbose = FALSE), "gaussian")
})

test_that("detect_family returns gaussian for integer vector with >2 unique values", {
    expect_equal(detect_family(c(0, 1, 2, 3), verbose = FALSE), "gaussian")
})

test_that("detect_family handles NAs correctly — ignores them in detection", {
    expect_equal(detect_family(c(0, 1, NA, 0, 1, NA), verbose = FALSE), "binomial")
})

test_that("detect_family logs detection result when verbose = TRUE", {
    msg <- capture.output(
        result <- detect_family(c(0, 1), verbose = TRUE),
        type = "message"
    )
    expect_equal(result, "binomial")
    # Verbose path shouldn't error even if logging goes to stderr
})

# ---------------------------------------------------------------------------
# recode_phenotype() tests
# ---------------------------------------------------------------------------

test_that("recode_phenotype converts {1,2} to {0,1}", {
    y_in  <- c(1, 2, 1, 2, 2)
    y_out <- recode_phenotype(y_in, verbose = FALSE)
    expect_equal(as.numeric(y_out), c(0, 1, 0, 1, 1))
})

test_that("recode_phenotype converts {-1,1} to {0,1}", {
    y_in  <- c(-1, 1, -1, 1, -1)
    y_out <- recode_phenotype(y_in, verbose = FALSE)
    expect_equal(as.numeric(y_out), c(0, 1, 0, 1, 0))
})

test_that("recode_phenotype leaves {0,1} unchanged", {
    y_in  <- c(0, 1, 1, 0, 0)
    y_out <- recode_phenotype(y_in, verbose = FALSE)
    expect_equal(as.numeric(y_out), y_in)
})

test_that("recode_phenotype preserves NAs in output", {
    y_in  <- c(1, 2, NA, 1, NA)
    y_out <- recode_phenotype(y_in, verbose = FALSE)
    expect_true(is.na(y_out[3]))
    expect_true(is.na(y_out[5]))
    expect_equal(y_out[1], 0)
    expect_equal(y_out[2], 1)
})

test_that("recode_phenotype stores original_coding attribute", {
    y_12  <- recode_phenotype(c(1, 2), verbose = FALSE)
    y_n11 <- recode_phenotype(c(-1, 1), verbose = FALSE)
    y_01  <- recode_phenotype(c(0, 1), verbose = FALSE)
    expect_equal(attr(y_12,  "original_coding"), "12")
    expect_equal(attr(y_n11, "original_coding"), "-11")
    expect_equal(attr(y_01,  "original_coding"), "01")
})

test_that("recode_phenotype errors on unrecognised coding", {
    expect_error(recode_phenotype(c(0, 1, 2), verbose = FALSE))
})

test_that("recode_phenotype respects explicit from_coding argument", {
    y_out <- recode_phenotype(c(1, 2), from_coding = "12", verbose = FALSE)
    expect_equal(as.numeric(y_out), c(0, 1))
})

# ---------------------------------------------------------------------------
# validate_binary_phenotype() tests
# ---------------------------------------------------------------------------

test_that("validate_binary_phenotype returns TRUE for valid binary phenotype", {
    phe <- data.frame(ID = 1:30, Y = c(rep(0, 15), rep(1, 15)))
    expect_true(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype returns FALSE when NAs are present", {
    phe <- data.frame(ID = 1:5, Y = c(0, 1, NA, 0, 1))
    expect_false(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype returns FALSE for continuous phenotype", {
    phe <- data.frame(ID = 1:5, Y = c(1.1, 2.2, 3.3, 4.4, 5.5))
    expect_false(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype returns FALSE for >2 unique values", {
    phe <- data.frame(ID = 1:4, Y = c(0, 1, 2, 3))
    expect_false(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype warns when too few cases", {
    phe <- data.frame(ID = 1:15, Y = c(rep(0, 12), rep(1, 3)))
    expect_warning(
        validate_binary_phenotype(phe, min_cases = 10, verbose = FALSE),
        regexp = "only 3 cases"
    )
})

test_that("validate_binary_phenotype warns when too few controls", {
    phe <- data.frame(ID = 1:15, Y = c(rep(0, 3), rep(1, 12)))
    expect_warning(
        validate_binary_phenotype(phe, min_controls = 10, verbose = FALSE),
        regexp = "only 3 controls"
    )
})

test_that("validate_binary_phenotype still returns TRUE with too-few warning", {
    phe <- data.frame(ID = 1:15, Y = c(rep(0, 12), rep(1, 3)))
    result <- suppressWarnings(
        validate_binary_phenotype(phe, min_cases = 10, verbose = FALSE)
    )
    expect_true(result)
})

test_that("validate_binary_phenotype works with {1,2} coding", {
    phe <- data.frame(ID = 1:20, Y = c(rep(1, 10), rep(2, 10)))
    expect_true(validate_binary_phenotype(phe, verbose = FALSE))
})
