#include "rMVP.h"

// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::mat GInv(const arma::mat& A){
	
	arma::mat ginv;
	if(A.n_rows == 1){
		ginv = 1 / A;
	}else{
		try {
			ginv = arma::inv_sympd(A);
		}
		catch (const std::exception& e) {
			arma::mat U;
			arma::vec s;
			arma::mat V;
			double tol = sqrt(datum::eps);
			
			svd(U,s,V,A);
			U = conv_to<mat>::from(conj(conv_to<cx_mat>::from(U)));
			arma::vec sMax(2); sMax.fill(0);
			sMax[1] = tol * s[0];
			arma::uvec Positive = find(s > sMax.max());
			arma::mat Up = U.cols(Positive);
			Up.each_row() %= 1/s(Positive).t();
			ginv = V.cols(Positive) * Up.t();
		}
	}
	return ginv;
}

template <typename T>
NumericVector getRow(const XPtr<BigMatrix> pMat, const int row){
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = pMat->ncol();

	NumericVector snp(ind);

	for(int i = 0; i < ind; i++){
		snp[i] = genomat[i][row];
	}

	return snp;
}

// [[Rcpp::export]]
NumericVector getRow(SEXP pBigMat, const int row){

	XPtr<BigMatrix> xpMat(pBigMat);

	DISPATCH_MATRIX_TYPE(getRow, xpMat, row);
}

template <typename T>
SEXP glm_c(const arma::vec &y, const arma::mat &X, const arma::mat & iXX, XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	
	int n;
	uvec _geno_ind;
	if(geno_ind.isNotNull()){
        _geno_ind = as<uvec>(geno_ind) - 1;
		n = _geno_ind.n_elem;
    }else{
		n = marker_bycol ? pMat->nrow() : pMat->ncol();
	}

	int m;
	uvec _marker_ind;
	if(marker_ind.isNotNull()){
        _marker_ind = as<uvec>(marker_ind) - 1;
		m = _marker_ind.n_elem;
    }else{
		m = marker_bycol ? pMat->ncol() : pMat->nrow();
	}

	int q0 = X.n_cols;
	if(y.n_elem != (arma::uword)n)	throw Rcpp::exception("number of individuals not match!");

	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

	arma::mat xy = X.t() * y;
	double yy = sum(y % y);
	arma::mat res(m, 1 + 1 + 1 + q0);

	arma::mat iXXs(q0 + 1, q0 + 1);

	arma::mat Z_buffer(n, step, fill::none);
	int i = 0, j = 0;
	int i_marker = 0;
	while(i < m) {
		
		int cnt = 0;
		for (; j < m && cnt < step; j++)
		{
			cnt++;
		}

		fill_geno_buffer(Z_buffer, genomat, cnt, i_marker, _geno_ind, _marker_ind, marker_bycol, true);

		#pragma omp parallel for firstprivate(iXXs)
		for(int l = 0; l < cnt; l++){
			arma::mat xs(q0, 1);
			arma::mat B21(1, q0);
			arma::mat NeginvB22B21(1, q0);
			arma::mat rhs(xy.n_rows + 1, 1);
			arma::mat beta(q0 + 1, 1);
			arma::vec se(q0 + 1);
			arma::vec pvalue(q0 + 1);

			double sy = sum(Z_buffer.col(l) % y);
			double ss = sum(Z_buffer.col(l) % Z_buffer.col(l));
			xs = X.t() * Z_buffer.col(l);
			B21 = xs.t() * iXX;
			double t2 = as_scalar(B21 * xs);
			double B22 = (ss - t2);
			double invB22;
			int df;
			if(B22 < 1e-8){
				invB22 = 0;
				df = n - q0;
			}else{
				invB22 = 1 / B22;
				df = n - q0 - 1;
			}
			NeginvB22B21 = -1 * invB22 * B21;
			iXXs(q0, q0) = invB22;
			iXXs.submat(0, 0, q0 - 1, q0 - 1) = iXX + invB22 * B21.t() * B21;
			iXXs(q0, arma::span(0, q0 - 1)) = NeginvB22B21;
			iXXs(arma::span(0, q0 - 1), q0) = NeginvB22B21.t();

			// statistics
			rhs.rows(0, xy.n_rows - 1) = xy;
			rhs(xy.n_rows, 0) = sy;
			beta = iXXs * rhs;

			double ve = (yy - as_scalar(beta.t() * rhs)) / df;
			for(int ff = 0; ff < (q0 + 1); ff++){
				se[ff] = sqrt(iXXs(ff, ff) * ve);
				pvalue[ff] = 2 * R::pt(abs(beta[ff] / se[ff]), df, false, false);
				res(i_marker + l, ff + 2) = pvalue[ff];
			}

			if(invB22 == 0){
				beta[q0] = NA_REAL;
				se[q0] = NA_REAL;
				res(i_marker + l, q0) = NA_REAL;
			}
			res(i_marker + l, 0) = beta[q0];
			res(i_marker + l, 1) = se[q0]; 
		}

		i = j;
		i_marker += cnt;
		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();
	return wrap(res);
}

// [[Rcpp::export]]
SEXP glm_c(const arma::vec & y, const arma::mat & X, const arma::mat & iXX, SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return glm_c<char>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 2:
		return glm_c<short>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 4:
		return glm_c<int>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 8:
		return glm_c<double>(y, X, iXX, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP mlm_c(const arma::vec & y, const arma::mat & X, const arma::mat & U, const double vgs, XPtr<BigMatrix> pMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	
	int n;
	uvec _geno_ind;
	if(geno_ind.isNotNull()){
        _geno_ind = as<uvec>(geno_ind) - 1;
		n = _geno_ind.n_elem;
    }else{
		n = marker_bycol ? pMat->nrow() : pMat->ncol();
	}

	int m;
	uvec _marker_ind;
	if(marker_ind.isNotNull()){
        _marker_ind = as<uvec>(marker_ind) - 1;
		m = _marker_ind.n_elem;
    }else{
		m = marker_bycol ? pMat->ncol() : pMat->nrow();
	}

	int q0 = X.n_cols;
	if(y.n_elem != (arma::uword)n)	throw Rcpp::exception("number of individuals not match.!");

	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

	arma::mat Uy = U.t() * y;
	arma::mat UX = U.t() * X;
	arma::mat UXUy = UX.t() * Uy;
	arma::mat iUXUX = GInv(UX.t() * UX);
	
	arma::mat res(m, 3);		
	arma::mat iXXs(q0 + 1, q0 + 1);

	arma::mat Z_buffer(n, step, fill::none);
	int i = 0, j = 0;
	int i_marker = 0;
	while(i < m) {
		
		int cnt = 0;
		for (; j < m && cnt < step; j++)
		{
			cnt++;
		}

		fill_geno_buffer(Z_buffer, genomat, cnt, i_marker, _geno_ind, _marker_ind, marker_bycol, true);

		#pragma omp parallel for firstprivate(iXXs)
		for(int l = 0; l < cnt; l++){
			arma::mat Us = U.t() * Z_buffer.col(l);
			arma::mat UXUs = UX.t() * Us;
			double UsUs = as_scalar(Us.t() * Us);
			double UsUy = as_scalar(Us.t() * Uy);
			double B22 = UsUs - as_scalar(UXUs.t() * iUXUX * UXUs);
			double invB22 = 1 / B22;
			arma::mat B21 = UXUs.t() * iUXUX;
			arma::mat NeginvB22B21 = -1 * invB22 * B21;
			
			iXXs(q0, q0) = invB22;
			iXXs.submat(0, 0, q0 - 1, q0 - 1) = iUXUX + invB22 * B21.t() * B21;
			iXXs(q0, arma::span(0, q0 - 1)) = NeginvB22B21;
			iXXs(arma::span(0, q0 - 1), q0) = NeginvB22B21.t();

			// statistics
			arma::mat rhs(UXUy.n_rows + 1, 1);
			rhs.rows(0, UXUy.n_rows - 1) = UXUy;
			rhs(UXUy.n_rows, 0) = UsUy;
			arma::mat beta = iXXs * rhs;
			int df = n - q0 - 1;

			res(i_marker + l, 0) = beta(q0, 0);
			res(i_marker + l, 1) = sqrt(iXXs(q0, q0) * vgs); 
			res(i_marker + l, 2) = 2 * R::pt(abs(res(i_marker + l, 0) / res(i_marker + l, 1)), df, false, false);
		}

		i = j;
		i_marker += cnt;
		if(!Progress::check_abort()) progress.increment(cnt);
	}
	Z_buffer.reset();
	return wrap(res);
}

// [[Rcpp::export]]
SEXP mlm_c(const arma::vec & y, const arma::mat & X, const arma::mat & U, const double vgs, SEXP pBigMat, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const bool marker_bycol = true, const int step = 10000, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return mlm_c<char>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 2:
		return mlm_c<short>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 4:
		return mlm_c<int>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	case 8:
		return mlm_c<double>(y, X, U, vgs, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

// ---------------------------------------------------------------------------
// Logistic regression C++ backend — Phase 2 (logistic-regression feature)
//
// Implements Newton-Raphson IRLS for each marker SNP with optional Firth
// penalised-likelihood correction for complete/quasi-separation stability.
// Wald test is used for significance (chi-squared with 1 df).
//
// Returns an m x 3 numeric matrix: [Effect, SE, p-value] per marker.
// Markers that fail to converge are returned as NA.
// ---------------------------------------------------------------------------

// Compute Firth correction score adjustment:
//   U_firth[j] = 0.5 * X[,j] * h_ii * (1 - 2*P_i)
// where h_ii = X[i,] (X'WX)^{-1} X[i,]'  (hat-matrix diagonal).
// Returns a (q+1)-length correction vector in coefficient space.
static arma::vec compute_firth_score(const arma::mat &X,
                                     const arma::vec &P,
                                     const arma::mat &iXtWX) {
    int n = X.n_rows;
    int p = X.n_cols;
    arma::vec h(n);
    for (int i = 0; i < n; i++) {
        arma::rowvec xi = X.row(i);
        h[i] = arma::as_scalar(xi * iXtWX * xi.t());
    }
    // Firth score: X' * diag(0.5*(1-2P)*h)
    arma::vec adj(p, arma::fill::zeros);
    for (int j = 0; j < p; j++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += X(i, j) * 0.5 * (1.0 - 2.0 * P[i]) * h[i];
        }
        adj[j] = s;
    }
    return adj;
}

// Core logistic_c template — called via SEXP dispatcher below.
template <typename T>
SEXP logistic_c(const arma::vec &y,
                const arma::mat &X,
                XPtr<BigMatrix> pMat,
                const Nullable<arma::uvec> geno_ind   = R_NilValue,
                const Nullable<arma::uvec> marker_ind = R_NilValue,
                const bool   marker_bycol = true,
                const int    step         = 10000,
                const bool   verbose      = true,
                const int    threads      = 0,
                const bool   firth        = true,
                const int    max_iter     = 25,
                const double tol          = 1e-8) {

    omp_setup(threads);

    MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

    // ---- resolve individual index ----------------------------------------
    int n;
    arma::uvec _geno_ind;
    if (geno_ind.isNotNull()) {
        _geno_ind = as<arma::uvec>(geno_ind) - 1;
        n = _geno_ind.n_elem;
    } else {
        n = marker_bycol ? pMat->nrow() : pMat->ncol();
    }

    // ---- resolve marker index --------------------------------------------
    int m;
    arma::uvec _marker_ind;
    if (marker_ind.isNotNull()) {
        _marker_ind = as<arma::uvec>(marker_ind) - 1;
        m = _marker_ind.n_elem;
    } else {
        m = marker_bycol ? pMat->ncol() : pMat->nrow();
    }

    int q0 = X.n_cols;   // number of null-model covariates (incl. intercept)
    if (y.n_elem != (arma::uword)n)
        throw Rcpp::exception("logistic_c: number of individuals does not match!");

    // ---- result storage: m rows x 3 cols (Effect, SE, p-value) ----------
    arma::mat res(m, 3, arma::fill::value(NA_REAL));

    MinimalProgressBar_plus pb;
    Progress progress(m, verbose, pb);

    // ---- genotype buffer -------------------------------------------------
    arma::mat Z_buffer(n, step, arma::fill::none);
    int i = 0, j = 0;
    int i_marker = 0;

    while (i < m) {
        // fill batch
        int cnt = 0;
        for (; j < m && cnt < step; j++) { cnt++; }
        fill_geno_buffer(Z_buffer, genomat, cnt, i_marker, _geno_ind, _marker_ind, marker_bycol, true);

        // ---- per-marker Newton-Raphson -----------------------------------
        #pragma omp parallel for schedule(dynamic)
        for (int l = 0; l < cnt; l++) {
            // Build augmented design matrix [X | z_l]
            arma::mat Xfull(n, q0 + 1);
            Xfull.cols(0, q0 - 1) = X;
            Xfull.col(q0)         = Z_buffer.col(l);

            // Check for monomorphic marker (zero variance) → skip
            double zvar = arma::var(Z_buffer.col(l));
            if (zvar < 1e-10) {
                // Leave as NA
                continue;
            }

            // Newton-Raphson IRLS
            arma::vec beta(q0 + 1, arma::fill::zeros);
            bool converged = false;

            for (int iter = 0; iter < max_iter; iter++) {
                arma::vec eta  = Xfull * beta;               // linear predictor
                arma::vec P    = 1.0 / (1.0 + arma::exp(-eta));  // P(Y=1)
                arma::vec W    = P % (1.0 - P);              // working weights

                // Clamp weights to avoid numerical zero (near-separation)
                W = arma::clamp(W, 1e-8, 0.25);

                // X' W X
                arma::mat XtW  = Xfull.t();
                XtW.each_row() %= W.t();
                arma::mat XtWX = XtW * Xfull;

                // Symmetric positive-definite solve (fallback to SVD on fail)
                arma::mat iXtWX;
                bool ok = arma::inv_sympd(iXtWX, XtWX);
                if (!ok) {
                    // SVD-based pseudo-inverse
                    arma::mat U2; arma::vec s; arma::mat V2;
                    arma::svd(U2, s, V2, XtWX);
                    double thr = arma::max(s) * 1e-10;
                    arma::vec si = 1.0 / arma::clamp(s, thr, arma::datum::inf);
                    iXtWX = V2 * arma::diagmat(si) * U2.t();
                }

                // Score: X'(Y - P) [+ Firth adjustment if enabled]
                arma::vec score = Xfull.t() * (y - P);
                if (firth) {
                    score += compute_firth_score(Xfull, P, iXtWX);
                }

                // Newton step: Δβ = (X'WX)^{-1} score
                arma::vec delta_beta = iXtWX * score;
                beta += delta_beta;

                if (arma::norm(delta_beta, "inf") < tol) {
                    converged = true;
                    break;
                }
            }

            if (!converged) {
                // Leave row as NA — marker failed convergence
                continue;
            }

            // ---- Final statistics ----------------------------------------
            arma::vec eta_f = Xfull * beta;
            arma::vec P_f   = 1.0 / (1.0 + arma::exp(-eta_f));
            arma::vec W_f   = arma::clamp(P_f % (1.0 - P_f), 1e-8, 0.25);

            arma::mat XtW_f  = Xfull.t();
            XtW_f.each_row() %= W_f.t();
            arma::mat XtWX_f = XtW_f * Xfull;

            arma::mat iXtWX_f;
            bool ok2 = arma::inv_sympd(iXtWX_f, XtWX_f);
            if (!ok2) {
                arma::mat U2; arma::vec s; arma::mat V2;
                arma::svd(U2, s, V2, XtWX_f);
                double thr = arma::max(s) * 1e-10;
                arma::vec si = 1.0 / arma::clamp(s, thr, arma::datum::inf);
                iXtWX_f = V2 * arma::diagmat(si) * U2.t();
            }

            double effect = beta[q0];
            double se     = std::sqrt(std::abs(iXtWX_f(q0, q0)));
            double z_stat = effect / se;
            // Wald chi-squared test: p = 2*Phi(-|z|)
            double pval   = 2.0 * R::pnorm(-std::abs(z_stat), 0.0, 1.0, 1, 0);

            res(i_marker + l, 0) = effect;
            res(i_marker + l, 1) = se;
            res(i_marker + l, 2) = pval;
        }

        i        = j;
        i_marker += cnt;
        if (!Progress::check_abort()) progress.increment(cnt);
    }

    Z_buffer.reset();
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
SEXP logistic_c(const arma::vec &y,
                const arma::mat &X,
                SEXP pBigMat,
                const Nullable<arma::uvec> geno_ind   = R_NilValue,
                const Nullable<arma::uvec> marker_ind = R_NilValue,
                const bool   marker_bycol = true,
                const int    step         = 10000,
                const bool   verbose      = true,
                const int    threads      = 0,
                const bool   firth        = true,
                const int    max_iter     = 25,
                const double tol          = 1e-8) {

    XPtr<BigMatrix> xpMat(pBigMat);

    switch (xpMat->matrix_type()) {
    case 1:
        return logistic_c<char>  (y, X, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads, firth, max_iter, tol);
    case 2:
        return logistic_c<short> (y, X, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads, firth, max_iter, tol);
    case 4:
        return logistic_c<int>   (y, X, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads, firth, max_iter, tol);
    case 8:
        return logistic_c<double>(y, X, xpMat, geno_ind, marker_ind, marker_bycol, step, verbose, threads, firth, max_iter, tol);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

