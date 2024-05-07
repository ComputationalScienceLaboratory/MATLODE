#include <algorithm>
#include <vector>

#include "cudaUtils.hpp"
#include "mex.h"

#pragma region GMRES Extra Algorithms
// TODO: replace with std::span
void prevgiv(double* column, const std::vector<double>& cs, const std::vector<double>& sn, int n_rot) {
	/* carry out n_rot rotations on column vector */

	int k;

	double h1, h2;

	for (k = 0; k < n_rot; k++) {
		h1            = column[k];
		h2            = column[k + 1];
		column[k]     = cs[k] * h1 + sn[k] * h2;
		column[k + 1] = -sn[k] * h1 + cs[k] * h2;
	}
	return;
}

void givens(double& cs, double& sn, double x, double y) {
	double zt;
	double z;
	double as;

	if (y == 0.0) {
		cs = 1.0;
		sn = 0.0;
	} else if (x == 0.0) {
		cs = 0.0;
		sn = 1.0;
	} else if (fabs(y) >= fabs(x)) {
		zt = x / y;
		z  = fabs(zt);
		as = 1.0 / sqrt(z * z + 1.0);
		cs = z * as;
		sn = as * (z / zt);
	} else {
		zt = y / x;
		z  = fabs(zt);
		cs = 1.0 / sqrt(z * z + 1);
		sn = cs * zt;
	}
}

void uptrisol(int neqns, std::vector<double>& hes_sol, const std::vector<double>& hes_mat, int col_length) {
	/* solve upper-triangular system */

	int k, j;

	for (k = neqns - 1; k >= 0; k--) {
		hes_sol[k] = hes_sol[k] / (hes_mat[k * col_length + k]);
		for (j = k - 1; j >= 0; j--) { hes_sol[j] -= hes_sol[k] * (hes_mat[k * col_length + j]); }
	}
	return;
}
#pragma endregion

DevPtr<double> gmresm(uint32_t m,                         /* restart frequency */
                      uint32_t maxnit,                    /* max number of iterations */
                      double tol,                         /* rel. convergence tolerance */
                      const SparseMatCSC<DevPtr>& matrix, /* matrix data structure */
                      DevPtr<double> rhs) {               /* right hand side vector */

	mxAssert(matrix.rows == matrix.cols, "Matrix must be square");
	mxAssert(matrix.rows == rhs.size, "Rhs vector dimension must match the matrix size");

	CudaSparseMath gpu_math;

	// Residual vector
	DevPtr<double> res = new_gpu_buffer<double>(rhs.size);
	// Solution vector (return value)
	DevPtr<double> sol = new_gpu_buffer<double>(rhs.size);

	cusparseSpMatDescr_t mat_cuda = cuda_csc(matrix);
	cusparseDnVecDescr_t rhs_cuda = cuda_vec(rhs);
	cusparseDnVecDescr_t sol_cuda = cuda_vec(sol);

	// Do setup work for matrix-vector products
	CHECK_CUSPARSE(gpu_math.prepare_sparse_matvec(mat_cuda, rhs_cuda, sol_cuda));

#pragma region GMRES Algorithm
	double rhs_nrm = CHECK_CUBLAS_VARIANT(gpu_math.norm(rhs));
	double res_tol = tol * rhs_nrm;

	// Scratch space
	DevPtr<double> tmp_vec        = new_gpu_buffer<double>(rhs.size);
	cusparseDnVecDescr_t tmp_cuda = cuda_vec(tmp_vec);

	// Calculate initial residual
	CHECK_CUBLAS(gpu_math.copy(res, rhs));

	CHECK_CUSPARSE(gpu_math.sparse_matvec(mat_cuda, sol_cuda, tmp_cuda));

	CHECK_CUBLAS(gpu_math.axpy(-1, tmp_vec, res));

	// Keep track of the size of the residual - this is what will let us know when we're ready to be done
	double res_nrm = CHECK_CUBLAS_VARIANT(gpu_math.norm(res));
	if (res_nrm <= res_tol) {
		// mexPrintf("CUDA GMRES exited early w/ residual norm %e < %e\n", res_nrm, res_tol);
		tmp_vec.free();
		res.free();

		CHECK_CUSPARSE(cusparseDestroySpMat(mat_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(rhs_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(sol_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(tmp_cuda));
		return sol;
	}

	int m1 = m + 1;
	// Space for hessenberg matrix
	std::vector<double> hes_mat(m1 * m), hes_rhs(m1), hes_sol(m), cs(m), sn(m);

	// For alignment reasons, rhs.size needs to be a multiple of 4 here, so we'll round up to the nearest multiple of 4
	size_t rhs_buffer_size = rhs.size + 4 - rhs.size % 4;
	DevPtr<double> vv      = new_gpu_buffer<double>(rhs_buffer_size * m1);

	uint32_t current_iter = 0;
	uint32_t inner_iters  = 0;
	while ((res_nrm > res_tol) && (current_iter < maxnit)) {
		current_iter += 1;
		gpu_math.copy(vv, res);
		double tmp = 1.0 / res_nrm;

		gpu_math.scale(vv, tmp);
		std::fill(hes_rhs.begin(), hes_rhs.end(), 0);
		hes_rhs[0] = res_nrm;

		size_t k;
		size_t k1;
		for (k = 0; k < m; k++) {
			k1 = k + 1;
			inner_iters += 1;

			DevPtr<double, true> k_vec_view  = vv.get_view(k * rhs_buffer_size, rhs.size),
			                     k1_vec_view = vv.get_view(k1 * rhs_buffer_size, rhs.size);

			cusparseDnVecDescr_t k_cuda  = cuda_vec(k_vec_view);
			cusparseDnVecDescr_t k1_cuda = cuda_vec(k1_vec_view);
			CHECK_CUSPARSE(gpu_math.sparse_matvec(mat_cuda, k_cuda, k1_cuda));
			CHECK_CUSPARSE(cusparseDestroyDnVec(k_cuda));
			CHECK_CUSPARSE(cusparseDestroyDnVec(k1_cuda));

			for (size_t j = 0; j < k1; j++) {
				DevPtr<double, true> j_vec_view = vv.get_view(j * rhs_buffer_size, rhs.size);

				hes_mat[k * m1 + j] = CHECK_CUBLAS_VARIANT(gpu_math.dot_product(k1_vec_view, j_vec_view));
				CHECK_CUBLAS(gpu_math.axpy(-hes_mat[k * m1 + j], j_vec_view, k1_vec_view));
			}
			hes_mat[k * m1 + k1] = CHECK_CUBLAS_VARIANT(gpu_math.norm(k1_vec_view));
			if (hes_mat[k * m1 + k1] == 0.0) { break; }
			tmp = 1.0 / hes_mat[k * m1 + k1];
			CHECK_CUBLAS(gpu_math.scale(k1_vec_view, tmp));

			prevgiv(&hes_mat[k * m1], cs, sn, k);
			double h1 = hes_mat[k * m1 + k];
			double h2 = hes_mat[k * m1 + k1];
			givens(cs[k], sn[k], h1, h2);
			hes_mat[k * m1 + k]  = cs[k] * h1 + sn[k] * h2;
			hes_mat[k * m1 + k1] = -sn[k] * h1 + cs[k] * h2;
			hes_rhs[k1]          = -sn[k] * hes_rhs[k];
			hes_rhs[k]           = cs[k] * hes_rhs[k];

			if (fabs(hes_rhs[k1]) <= res_tol) break;
		}

		size_t meff = k1;
		std::copy(hes_rhs.begin(), hes_rhs.end() - 1, hes_sol.begin());
		uptrisol(meff, hes_sol, hes_mat, m1);

		for (k = 0; k < meff; k++) {
			CHECK_CUBLAS(gpu_math.axpy(hes_sol[k], vv.get_view(k * rhs_buffer_size, rhs.size), sol));
		}
		CHECK_CUBLAS(gpu_math.copy(res, rhs));
		CHECK_CUSPARSE(gpu_math.sparse_matvec(mat_cuda, sol_cuda, tmp_cuda));
		CHECK_CUBLAS(gpu_math.axpy(-1.0, tmp_vec, res));

		res_nrm = CHECK_CUBLAS_VARIANT(gpu_math.norm(res));
	}
#pragma endregion

	// mexPrintf("cudaGmres converged in %d inner iterations with a residual norm of %e\n", inner_iters, res_nrm);

	vv.free();
	tmp_vec.free();
	res.free();

	CHECK_CUSPARSE(cusparseDestroySpMat(mat_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(rhs_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(sol_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(tmp_cuda));

	return sol;
}

std::optional<uint32_t> get_uint(const mxArray* input) {
	if (mxIsDouble(input)) {
		double in = mxGetDoubles(input)[0];
		if (std::trunc(in) == in) {
			return in;
		} else {
			return {};
		}
	} else if (mxIsInt16(input)) {
		return mxGetInt16s(input)[0];
	} else if (mxIsInt32(input)) {
		return mxGetInt32s(input)[0];
	} else if (mxIsUint16(input)) {
		return mxGetUint16s(input)[0];
	} else if (mxIsUint32(input)) {
		return mxGetUint32s(input)[0];
	} else {
		return {};
	}
}

void checkArguments(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	if (nrhs < 2) { mexErrMsgIdAndTxt("Matlode:cudaGmres:nrhs", "Two inputs required"); }

	if (nlhs != 1) { mexErrMsgIdAndTxt("Matlode:cudaGmres:nlhs", "One output required"); }

	// First Argument - A (m-by-n matrix of doubles)
	// TODO: make work with floats as well

	const size_t* a_dimensions = mxGetDimensions(prhs[0]);
	size_t a_num_dimension     = mxGetNumberOfDimensions(prhs[0]);
	if (a_num_dimension != 2) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:a_dims", "Input 1 (A) must be n-by-n dimension (matrix)");
	}

	if (a_dimensions[0] != a_dimensions[1]) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:a_dims", "Input 1 (A) must be square matrix");
	}

	if (!mxIsSparse(prhs[0])) { mexErrMsgIdAndTxt("Matlode:cudaGmres:a_sparse", "Input 1 (A) must be sparse"); }

	// Second Argument - b (n-row vector of doubles)

	if (mxGetNumberOfDimensions(prhs[1]) > 2 ||
	    mxGetNumberOfDimensions(prhs[1]) == 2 &&
	        (mxGetDimensions(prhs[1])[0] != 1 && mxGetDimensions(prhs[1])[1] != 1)) {
		std::ostringstream ss;
		ss << "Input 2 (b) must be a vector. Instead found ";
		size_t i;
		for (i = 0; i < mxGetNumberOfDimensions(prhs[1]) - 1; i++) { ss << mxGetDimensions(prhs[1])[i] << 'x'; }
		ss << mxGetDimensions(prhs[1])[i];
		mexErrMsgIdAndTxt("Matlode:cudaGmres:b_dim", ss.str().c_str());
	}

	if (mxGetNumberOfElements(prhs[1]) != a_dimensions[1]) {
		std::ostringstream ss;
		ss << "Input 2 (b) length must match Input 1 (A) matrix dimension. " << mxGetNumberOfElements(prhs[1])
		   << " v.s. " << a_dimensions[0] << 'x' << a_dimensions[1];
		mexErrMsgIdAndTxt("Matlode:cudaGmres:b_dim", ss.str().c_str());
	}

	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:b_type", "Input 2 (b) must be noncomplex double");
	}

	// Optional arguments

	if (nrhs < 3) { return; }
	// Third argument - restart

	if (mxGetNumberOfElements(prhs[2]) != 1) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:restart_dim", "Input 3 (restart) must be a scalar");
	}

	if (!get_uint(prhs[2]).has_value()) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:restart_type", "Input 3 (restart) must be convertible to Uint32");
	}

	if (nrhs < 4) { return; }
	// Fourth argument - tolerance

	if (mxGetNumberOfElements(prhs[3]) != 1) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:tol_dim", "Input 4 (tolerance) must be a scalar");
	}

	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:tol_type", "Input 4 (tolerance) must be noncomplex double");
	}

	if (mxGetDoubles(prhs[3])[0] <= 0) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:tol_pos", "Input 4 (tolerance) must be positive");
	}

	if (nrhs < 5) { return; }
	// Fifth argument - maximum iterations

	if (mxGetNumberOfElements(prhs[4]) != 1) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:maxit_dim", "Input 5 (maxit) must be a scalar");
	}

	if (!get_uint(prhs[4]).has_value()) {
		mexErrMsgIdAndTxt("Matlode:cudaGmres:maxit_type", "Input 5 (maxit) must be convertible to Uint32");
	}
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	checkArguments(nlhs, plhs, nrhs, prhs);

	// matrix A
	SparseMatCSC<HostPtr> a_host = get_csc(prhs[0]);
	SparseMatCSC<DevPtr> a_dev   = move_csc(a_host);

	// Vector b
	HostPtr<double> b_host = {mxGetDoubles(prhs[1]), a_host.cols, true};
	DevPtr<double> b_dev   = new_gpu_buffer(b_host);

	uint32_t restart = 0;
	if (nrhs >= 3) { restart = get_uint(prhs[2]).value_or(restart); }

	double tol = 1e-6;
	if (nrhs >= 4) { tol = mxGetDoubles(prhs[3])[0]; }

	uint32_t maxit = std::max(a_host.cols > UINT32_MAX ? UINT32_MAX : (uint32_t) a_host.cols, 10u);
	if (nrhs >= 5) { maxit = get_uint(prhs[4]).value_or(maxit); }

	if (restart == 0) { restart = maxit; }

	DevPtr<double> sol_dev = gmresm(restart, maxit, tol, a_dev, b_dev);

	// Output vector
	plhs[0]             = mxCreateDoubleMatrix(a_host.rows, 1, mxREAL);
	HostPtr<double> out = {mxGetDoubles(plhs[0]), a_host.rows, true};

	CHECK_CUDA(cuda_memcpy(out, sol_dev));

	b_dev.free();
	sol_dev.free();
}