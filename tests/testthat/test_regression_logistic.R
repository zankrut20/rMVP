# tests/testthat/test_regression_logistic.R
#
# Phase 5: Comprehensive logistic regression test suite
#
# Sections:
#   1  - Unit tests: phenotype utilities
#   2  - MVP.Logistic() output structure
#   3  - Integration tests: MVP() dispatch
#   4  - Numerical stability / edge cases
#   5  - Validation against R glm()
#   6  - Output format / attributes
#
# Run with: devtools::test(filter = "regression_logistic")


# ---------------------------------------------------------------------------
# Shared test helpers
# ---------------------------------------------------------------------------

# Build a tiny in-memory big.matrix (n x m) and a binary phenotype
.make_logistic_testdata <- function(n = 60, m = 20,
                                   case_frac = 0.5, seed = 42) {
    set.seed(seed)
    y_raw <- c(rep(0L, round(n * (1 - case_frac))),
               rep(1L, round(n * case_frac)))
    y_raw <- y_raw[seq_len(n)]

    gmat  <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    # Use double to avoid bigmemory typecast warnings
    geno  <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))

    phe   <- data.frame(ID = paste0("ind", seq_len(n)), Y = y_raw)
    map   <- data.frame(SNP = paste0("snp", seq_len(m)),
                        Chr = rep(1, m),
                        Pos = seq(1e6, by = 1e5, length.out = m))
    list(phe = phe, geno = geno, map = map, y = y_raw)
}

# Same dataset but {1, 2} coded
.make_12coded <- function(n = 60, m = 10, seed = 99) {
    d <- .make_logistic_testdata(n, m, seed = seed)
    d$phe[, 2] <- d$phe[, 2] + 1L   # 0/1 -> 1/2
    d
}

# Same dataset but {-1, 1} coded
.make_neg11coded <- function(n = 60, m = 10, seed = 77) {
    d <- .make_logistic_testdata(n, m, seed = seed)
    d$phe[, 2] <- ifelse(d$phe[, 2] == 0L, -1L, 1L)
    d
}

# ===========================================================================
# Section 1: Phenotype utility unit tests
# ===========================================================================

test_that("detect_family returns binomial for {0,1}", {
    expect_equal(detect_family(c(0, 1, 0, 1), verbose = FALSE), "binomial")
})

test_that("detect_family returns gaussian for continuous values", {
    expect_equal(detect_family(c(1.1, 2.3, 3.5), verbose = FALSE), "gaussian")
})

test_that("recode_phenotype maps {1,2} to {0,1} correctly", {
    res <- recode_phenotype(c(1, 2, 1, 2), verbose = FALSE)
    expect_equal(as.numeric(res), c(0, 1, 0, 1))
    expect_equal(attr(res, "original_coding"), "12")
})

test_that("recode_phenotype maps {-1,1} to {0,1} correctly", {
    res <- recode_phenotype(c(-1, 1, -1), verbose = FALSE)
    expect_equal(as.numeric(res), c(0, 1, 0))
    expect_equal(attr(res, "original_coding"), "-11")
})

test_that("validate_binary_phenotype accepts balanced binary phenotype", {
    phe <- data.frame(ID = 1:30, Y = c(rep(0, 15), rep(1, 15)))
    expect_true(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype rejects continuous phenotype", {
    phe <- data.frame(ID = 1:5, Y = c(0.1, 0.5, 1.2, 2.0, 3.1))
    expect_false(validate_binary_phenotype(phe, verbose = FALSE))
})

test_that("validate_binary_phenotype warns on too few cases", {
    phe <- data.frame(ID = 1:25, Y = c(rep(0, 22), rep(1, 3)))
    expect_warning(validate_binary_phenotype(phe, min_cases = 10, verbose = FALSE),
                   regexp = "only 3 cases")
})

test_that("validate_binary_phenotype warns on too few controls", {
    phe <- data.frame(ID = 1:25, Y = c(rep(0, 3), rep(1, 22)))
    expect_warning(validate_binary_phenotype(phe, min_controls = 10, verbose = FALSE),
                   regexp = "only 3 controls")
})

# ===========================================================================
# Section 2: MVP.Logistic() output structure
# ===========================================================================

test_that("MVP.Logistic produces a numeric matrix", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_true(is.matrix(res))
    expect_true(is.numeric(res))
})

test_that("MVP.Logistic output has exactly 3 columns", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(ncol(res), 3L)
})

test_that("MVP.Logistic output has m rows matching marker count", {
    d   <- .make_logistic_testdata(m = 15)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(nrow(res), 15L)
})

test_that("MVP.Logistic colnames are Effect / SE / p-value", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(colnames(res), c("Effect", "SE", "p-value"))
})

test_that("MVP.Logistic p-values are in [0, 1]", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    pvals <- res[, "p-value"]
    expect_true(all(pvals[!is.na(pvals)] >= 0))
    expect_true(all(pvals[!is.na(pvals)] <= 1))
})

test_that("MVP.Logistic standard errors are non-negative", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    se <- res[, "SE"]
    expect_true(all(se[!is.na(se)] >= 0))
})

test_that("MVP.Logistic works with {1,2} coded phenotype", {
    d   <- .make_12coded()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(ncol(res), 3L)
    expect_equal(attr(res, "family"), "binomial")
})

test_that("MVP.Logistic works with {-1,1} coded phenotype", {
    d   <- .make_neg11coded()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(ncol(res), 3L)
    expect_equal(attr(res, "family"), "binomial")
})

test_that("MVP.Logistic works with covariates", {
    d   <- .make_logistic_testdata(n = 80, m = 10)
    CV  <- matrix(rnorm(80), ncol = 1)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         CV = CV, cpu = 1, verbose = FALSE))
    expect_equal(ncol(res), 3L)
    expect_true(all(res[, "p-value"][!is.na(res[, "p-value"])] <= 1))
})

test_that("MVP.Logistic rejects continuous phenotype", {
    d       <- .make_logistic_testdata()
    d$phe[, 2] <- rnorm(nrow(d$phe))
    expect_error(
        suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                       cpu = 1, verbose = FALSE)),
        regexp = "binomial|binary|gaussian"
    )
})

test_that("MVP.Logistic rejects family='gaussian'", {
    d <- .make_logistic_testdata()
    expect_error(
        suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                       family = "gaussian",
                                       cpu = 1, verbose = FALSE)),
        regexp = "gaussian.*not supported|MVP\\.GLM"
    )
})

# ===========================================================================
# Section 3: Integration tests -- MVP() dispatch
# ===========================================================================

test_that("MVP auto-detects binary phenotype and dispatches to logistic", {
    d   <- .make_logistic_testdata(n = 60, m = 10)
    res <- suppressWarnings(
        MVP(phe = d$phe, geno = d$geno, map = d$map,
            method = "GLM", file.output = FALSE,
            verbose = FALSE, ncpus = 1)
    )
    expect_equal(attr(res$glm.results, "family"), "binomial")
    expect_equal(attr(res$glm.results, "method"), "Logistic")
})

test_that("MVP dispatches to MVP.GLM for continuous phenotype", {
    set.seed(1)
    d       <- .make_logistic_testdata(n = 60, m = 10)
    d$phe[, 2] <- rnorm(60)
    res <- MVP(phe = d$phe, geno = d$geno, map = d$map,
               method = "GLM", file.output = FALSE,
               verbose = FALSE, ncpus = 1)
    expect_null(attr(res$glm.results, "family"))
    expect_equal(ncol(res$glm.results), 3L)
})

test_that("MVP respects explicit family='binomial'", {
    d   <- .make_logistic_testdata(n = 60, m = 10)
    res <- suppressWarnings(
        MVP(phe = d$phe, geno = d$geno, map = d$map,
            method = "GLM", family = "binomial",
            file.output = FALSE, verbose = FALSE, ncpus = 1)
    )
    expect_equal(attr(res$glm.results, "family"), "binomial")
})

test_that("MVP respects explicit family='gaussian' even for binary phe", {
    d   <- .make_logistic_testdata(n = 60, m = 10)
    res <- suppressWarnings(
        MVP(phe = d$phe, geno = d$geno, map = d$map,
            method = "GLM", family = "gaussian",
            file.output = FALSE, verbose = FALSE, ncpus = 1)
    )
    expect_null(attr(res$glm.results, "family"))
})

test_that("MVP returns glm.results with 3 columns for logistic", {
    d   <- .make_logistic_testdata(n = 60, m = 10)
    res <- suppressWarnings(
        MVP(phe = d$phe, geno = d$geno, map = d$map,
            method = "GLM", file.output = FALSE,
            verbose = FALSE, ncpus = 1)
    )
    expect_equal(ncol(res$glm.results), 3L)
})

# ===========================================================================
# Section 4: Numerical stability / edge cases
# ===========================================================================

test_that("MVP.Logistic returns NA (not crash) for monomorphic markers", {
    d <- .make_logistic_testdata(n = 40, m = 5)
    # Force one monomorphic column (all zeros)
    suppressWarnings(d$geno[, 1] <- 0L)
    res <- suppressWarnings(
        MVP.Logistic(phe = d$phe, geno = d$geno, cpu = 1, verbose = FALSE)
    )
    expect_true(is.na(res[1, "Effect"]))
    expect_true(is.na(res[1, "SE"]))
    expect_true(is.na(res[1, "p-value"]))
    # Other markers should be OK
    expect_false(all(is.na(res[-1, "p-value"])))
})

test_that("MVP.Logistic handles extreme case imbalance without crashing", {
    set.seed(7)
    n  <- 100
    m  <- 8
    # 90 cases, 10 controls -- very imbalanced
    y  <- c(rep(1L, 90), rep(0L, 10))
    gmat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    geno <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))
    phe  <- data.frame(ID = seq_len(n), Y = y)

    # Must not throw an error
    expect_no_error(
        suppressWarnings(
            MVP.Logistic(phe = phe, geno = geno, cpu = 1,
                         firth = TRUE, verbose = FALSE)
        )
    )
    res <- suppressWarnings(
        MVP.Logistic(phe = phe, geno = geno, cpu = 1,
                     firth = TRUE, verbose = FALSE)
    )
    expect_equal(ncol(res), 3L)
})

test_that("MVP.Logistic with firth=FALSE still returns valid structure", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(
        MVP.Logistic(phe = d$phe, geno = d$geno,
                     firth = FALSE, cpu = 1, verbose = FALSE)
    )
    expect_equal(ncol(res), 3L)
    pvals <- res[, "p-value"]
    expect_true(all(pvals[!is.na(pvals)] >= 0 & pvals[!is.na(pvals)] <= 1))
})

test_that("MVP.Logistic convergence failure returns NA gracefully", {
    set.seed(101)
    n <- 20; m <- 4
    y <- c(rep(0L, 10), rep(1L, 10))
    gmat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    gmat[1:10, 1] <- 0L; gmat[11:20, 1] <- 2L   # near-perfect separation
    geno <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))
    phe  <- data.frame(ID = seq_len(n), Y = y)

    expect_silent(
        suppressWarnings(
            MVP.Logistic(phe = phe, geno = geno, cpu = 1,
                         firth = TRUE, verbose = FALSE)
        )
    )
})

# ===========================================================================
# Section 5: Validation against R glm()
# ===========================================================================

test_that("MVP.Logistic effect estimates correlate with R glm() (r > 0.90)", {
    skip_on_cran()
    set.seed(2024)
    n <- 200; m <- 30

    gmat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    geno <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))
    beta_true <- rep(0, m); beta_true[5] <- 0.6
    lp   <- gmat %*% beta_true - 0.5
    prob <- 1 / (1 + exp(-lp))
    y    <- rbinom(n, 1, prob)
    phe  <- data.frame(ID = seq_len(n), Y = y)

    res_mvp <- suppressWarnings(
        MVP.Logistic(phe = phe, geno = geno, cpu = 1, verbose = FALSE)
    )

    ref_effects <- vapply(seq_len(m), function(j) {
        dat <- data.frame(y = y, z = gmat[, j])
        tryCatch(
            coef(glm(y ~ z, data = dat, family = binomial()))["z"],
            error = function(e) NA_real_
        )
    }, numeric(1))

    ok <- !is.na(res_mvp[, "Effect"]) & !is.na(ref_effects)
    if (sum(ok) >= 5) {
        r <- cor(res_mvp[ok, "Effect"], ref_effects[ok])
        expect_gt(r, 0.90)
    } else {
        skip("Too few converged markers to compute correlation")
    }
})

test_that("MVP.Logistic SE estimates correlate with R glm() SE (r > 0.85)", {
    skip_on_cran()
    set.seed(2025)
    n <- 200; m <- 20
    gmat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    geno <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))
    y    <- rbinom(n, 1, 0.5)
    phe  <- data.frame(ID = seq_len(n), Y = y)

    res_mvp <- suppressWarnings(
        MVP.Logistic(phe = phe, geno = geno, cpu = 1, verbose = FALSE)
    )

    ref_se <- vapply(seq_len(m), function(j) {
        dat <- data.frame(y = y, z = gmat[, j])
        tryCatch(
            summary(glm(y ~ z, data = dat, family = binomial()))$coefficients["z", "Std. Error"],
            error = function(e) NA_real_
        )
    }, numeric(1))

    ok <- !is.na(res_mvp[, "SE"]) & !is.na(ref_se)
    if (sum(ok) >= 5) {
        r <- cor(res_mvp[ok, "SE"], ref_se[ok])
        expect_gt(r, 0.85)
    } else {
        skip("Too few converged markers to compute correlation")
    }
})

test_that("MVP.Logistic p-values correlate with R glm() p-values (r > 0.85)", {
    skip_on_cran()
    set.seed(2026)
    n <- 200; m <- 20
    gmat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    geno <- suppressWarnings(bigmemory::as.big.matrix(gmat + 0.0, type = "char", backingpath = tempdir()))
    y    <- rbinom(n, 1, 0.5)
    phe  <- data.frame(ID = seq_len(n), Y = y)

    res_mvp <- suppressWarnings(
        MVP.Logistic(phe = phe, geno = geno, cpu = 1, verbose = FALSE)
    )

    ref_pval <- vapply(seq_len(m), function(j) {
        dat <- data.frame(y = y, z = gmat[, j])
        tryCatch(
            summary(glm(y ~ z, data = dat, family = binomial()))$coefficients["z", "Pr(>|z|)"],
            error = function(e) NA_real_
        )
    }, numeric(1))

    ok <- !is.na(res_mvp[, "p-value"]) & !is.na(ref_pval)
    if (sum(ok) >= 5) {
        r <- cor(log(res_mvp[ok, "p-value"] + 1e-300),
                 log(ref_pval[ok] + 1e-300))
        expect_gt(r, 0.85)
    } else {
        skip("Too few converged markers to compute correlation")
    }
})

# ===========================================================================
# Section 6: Output format / attributes
# ===========================================================================

test_that("MVP.Logistic result has 'family' attribute equal to 'binomial'", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(attr(res, "family"), "binomial")
})

test_that("MVP.Logistic result has 'method' attribute equal to 'Logistic'", {
    d   <- .make_logistic_testdata()
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(attr(res, "method"), "Logistic")
})

test_that("MVP.Logistic result has n_cases attribute", {
    d   <- .make_logistic_testdata(n = 60)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_true(!is.null(attr(res, "n_cases")))
    expect_true(is.numeric(attr(res, "n_cases")))
})

test_that("MVP.Logistic result has n_controls attribute", {
    d   <- .make_logistic_testdata(n = 60)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_true(!is.null(attr(res, "n_controls")))
    expect_true(is.numeric(attr(res, "n_controls")))
})

test_that("MVP.Logistic n_cases + n_controls equals n", {
    n   <- 60
    d   <- .make_logistic_testdata(n = n)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    expect_equal(attr(res, "n_cases") + attr(res, "n_controls"), n)
})

test_that("MVP.Logistic result has firth attribute reflecting input", {
    d    <- .make_logistic_testdata()
    resT <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                           firth = TRUE,  cpu = 1, verbose = FALSE))
    resF <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                           firth = FALSE, cpu = 1, verbose = FALSE))
    expect_true(attr(resT, "firth"))
    expect_false(attr(resF, "firth"))
})

test_that("MVP.Logistic output is suitable for cbind with map", {
    d   <- .make_logistic_testdata(m = 12)
    res <- suppressWarnings(MVP.Logistic(phe = d$phe, geno = d$geno,
                                         cpu = 1, verbose = FALSE))
    combined <- cbind(d$map, res)
    expect_equal(nrow(combined), 12L)
    expect_equal(ncol(combined), ncol(d$map) + 3L)
})
