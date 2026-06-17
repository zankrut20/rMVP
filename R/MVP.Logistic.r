# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Genome-wide Association Study via Logistic Regression
#'
#' Performs GWAS for binary phenotypes using logistic regression with
#' Newton-Raphson (IRLS) estimation and optional Firth penalised-likelihood
#' correction for numerical stability under complete or quasi-separation.
#'
#' Build date: 2026-06-17
#' Last update: 2026-06-17
#'
#' @author Zankrut Goyani (Phase 3, logistic-regression feature)
#'
#' @param phe Phenotype: an \eqn{n \times 2}{n x 2} matrix or data frame with
#'   columns \code{[ID, Y]}.  Supported binary codings: \{0,1\}, \{1,2\},
#'   \{-1,1\}.  When \code{family = "auto"} the coding is detected
#'   automatically; continuous phenotypes are rejected with an informative
#'   error.
#' @param geno Genotype \code{big.matrix} object (memory-mapped).  Can be
#'   stored as n×m (individuals by markers) when \code{mrk_bycol = TRUE} or
#'   m×n otherwise.
#' @param CV Optional covariate matrix (\eqn{n \times k}{n x k}), added as
#'   fixed effects alongside the intercept in the null model.
#' @param ind_idx Integer vector. Row/column indices of individuals to include.
#'   \code{NULL} (default) uses all individuals.
#' @param mrk_idx Integer vector. Row/column indices of markers to scan.
#'   \code{NULL} (default) scans all markers.
#' @param mrk_bycol Logical. \code{TRUE} (default) means markers are stored
#'   as columns (n×m layout).
#' @param maxLine Integer. Number of markers loaded into RAM at once. Smaller
#'   values reduce peak memory use. Default \code{5000}.
#' @param cpu Integer. Number of threads for parallel marker scanning.
#'   Default \code{1}.
#' @param family Character. Phenotype family: \code{"auto"} (default) detects
#'   binary vs. continuous; \code{"binomial"} forces logistic regression.
#'   Passing \code{"gaussian"} here is an error — use \code{MVP.GLM()} instead.
#' @param firth Logical. Apply Firth penalised-likelihood correction (default
#'   \code{TRUE}).  Recommended for studies with rare variants or small sample
#'   sizes.
#' @param max_iter Integer. Maximum Newton-Raphson iterations per marker.
#'   Default \code{25}.
#' @param tol Numeric. Convergence tolerance (infinity norm of coefficient
#'   update). Default \code{1e-8}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \eqn{m \times 3}{m x 3} numeric matrix with columns:
#'   \describe{
#'     \item{Effect}{Log-odds effect estimate (log-OR) at the tested allele.}
#'     \item{SE}{Standard error of the effect estimate.}
#'     \item{p-value}{Wald test p-value (two-sided).}
#'   }
#'   Markers that fail to converge or are monomorphic are returned as \code{NA}.
#'   The matrix carries the following attributes:
#'   \describe{
#'     \item{family}{\code{"binomial"}}
#'     \item{method}{\code{"Logistic"}}
#'     \item{n_cases}{Number of case individuals (Y = 1 after recoding).}
#'     \item{n_controls}{Number of control individuals (Y = 0 after recoding).}
#'     \item{firth}{Logical; whether Firth correction was applied.}
#'   }
#'
#' @export
#'
#' @examples
#' \donttest{
#' phePath  <- system.file("extdata", "07_other", "mvp.phe", package = "rMVP")
#' phenotype <- read.table(phePath, header = TRUE)
#' # Simulate a binary phenotype from the first trait column
#' phenotype[, 2] <- ifelse(phenotype[, 2] > median(phenotype[, 2], na.rm=TRUE), 1, 0)
#' idx <- !is.na(phenotype[, 2])
#' phenotype <- phenotype[idx, ]
#'
#' genoPath <- system.file("extdata", "06_mvp-impute", "mvp.imp.geno.desc", package = "rMVP")
#' genotype <- attach.big.matrix(genoPath)
#' genotype <- deepcopy(genotype, rows = idx)
#'
#' result <- MVP.Logistic(phe = phenotype, geno = genotype, cpu = 1)
#' str(result)
#' }
MVP.Logistic <- function(
    phe,
    geno,
    CV         = NULL,
    ind_idx    = NULL,
    mrk_idx    = NULL,
    mrk_bycol  = TRUE,
    maxLine    = 5000,
    cpu        = 1,
    family     = "auto",
    firth      = TRUE,
    max_iter   = 25,
    tol        = 1e-8,
    verbose    = TRUE
) {

    # ------------------------------------------------------------------
    # 1. Input validation
    # ------------------------------------------------------------------
    if (!is.big.matrix(geno))
        stop("MVP.Logistic: 'geno' must be a big.matrix object.")

    if (ncol(phe) < 2)
        stop("MVP.Logistic: 'phe' must have at least 2 columns [ID, Y].")

    # Resolve sample size
    if (is.null(ind_idx)) {
        n <- if (mrk_bycol) nrow(geno) else ncol(geno)
        if (nrow(phe) != n)
            stop("MVP.Logistic: number of individuals differs between 'phe' and 'geno'.")
    } else {
        n <- length(ind_idx)
        if (nrow(phe) != n)
            stop("MVP.Logistic: length of 'ind_idx' differs from nrow(phe).")
    }

    if (!is.null(CV) && nrow(CV) != n)
        stop("MVP.Logistic: nrow(CV) must equal the number of individuals.")

    # ------------------------------------------------------------------
    # 2. Phenotype processing
    # ------------------------------------------------------------------
    y_raw <- as.numeric(phe[, 2])

    # Reject explicit gaussian request
    if (identical(family, "gaussian"))
        stop("MVP.Logistic: family='gaussian' is not supported. Use MVP.GLM() for continuous phenotypes.")

    # Auto-detect or validate family
    detected_family <- detect_family(y_raw, verbose = verbose)
    if (identical(family, "auto")) {
        family <- detected_family
    }
    if (!identical(family, "binomial")) {
        stop("MVP.Logistic: phenotype does not appear to be binary. ",
             "Detected family: '", detected_family, "'. ",
             "Use MVP.GLM() for continuous phenotypes.")
    }

    # Validate binary phenotype (warns if too few cases/controls)
    if (!validate_binary_phenotype(phe, verbose = verbose))
        stop("MVP.Logistic: binary phenotype validation failed (see messages above).")

    # Recode to {0, 1}
    y <- recode_phenotype(y_raw, verbose = verbose)
    n_cases    <- sum(y == 1L, na.rm = TRUE)
    n_controls <- sum(y == 0L, na.rm = TRUE)

    # After recoding, NAs must not be present (validate_binary_phenotype already checks)
    if (anyNA(y))
        stop("MVP.Logistic: NA values present in phenotype after recoding.")

    # ------------------------------------------------------------------
    # 3. Design matrix (null model: intercept + covariates)
    # ------------------------------------------------------------------
    if (is.null(CV)) {
        X0 <- matrix(1.0, nrow = n, ncol = 1)
    } else {
        if (anyNA(CV))
            stop("MVP.Logistic: NA values are not allowed in covariates (CV).")
        # Drop constant columns
        cv_vary <- apply(CV, 2, function(x) length(unique(x)) > 1)
        CV <- CV[, cv_vary, drop = FALSE]
        X0 <- cbind(matrix(1.0, n, 1), CV)
    }

    # ------------------------------------------------------------------
    # 4. Call C++ backend
    # ------------------------------------------------------------------
    logging.log("Logistic regression scanning...\n", verbose = verbose)

    mkl_env({
        results <- logistic_c(
            y            = y,
            X            = X0,
            pBigMat      = geno@address,
            geno_ind     = ind_idx,
            marker_ind   = mrk_idx,
            marker_bycol = mrk_bycol,
            step         = maxLine,
            verbose      = verbose,
            threads      = cpu,
            firth        = firth,
            max_iter     = max_iter,
            tol          = tol
        )
    })

    # ------------------------------------------------------------------
    # 5. Format output  (m x 3: Effect, SE, p-value)
    # ------------------------------------------------------------------
    results <- as.matrix(results)
    colnames(results) <- c("Effect", "SE", "p-value")

    # Attach metadata attributes
    attr(results, "family")     <- "binomial"
    attr(results, "method")     <- "Logistic"
    attr(results, "n_cases")    <- n_cases
    attr(results, "n_controls") <- n_controls
    attr(results, "firth")      <- firth

    return(results)
}
