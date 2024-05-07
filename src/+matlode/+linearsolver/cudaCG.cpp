#include <algorithm>
#include <vector>

#include "cudaUtils.hpp"
#include "mex.h"

DevPtr<double> conj_grad(uint32_t maxnit,                    /* max number of iterations */
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

#pragma region Conjugate Gradient Algorithm
	double rhs_nrm = CHECK_CUBLAS_VARIANT(gpu_math.norm(rhs));
	double res_tol = tol * rhs_nrm;

	// Scratch space
	DevPtr<double> tmp_vec        = new_gpu_buffer<double>(rhs.size);
	cusparseDnVecDescr_t tmp_cuda = cuda_vec(tmp_vec);

	// Keep track of the previous residual
	DevPtr<double> prev_res = new_gpu_buffer<double>(res.size);

	// Calculate initial residual
	CHECK_CUBLAS(gpu_math.copy(res, rhs));
	CHECK_CUSPARSE(gpu_math.sparse_matvec(mat_cuda, sol_cuda, tmp_cuda));
	CHECK_CUBLAS(gpu_math.axpy(-1, tmp_vec, res));

	// Keep track of the size of the residual - this is what will let us know when we're ready to be done
	double res_nrm = CHECK_CUBLAS_VARIANT(gpu_math.norm(res));
	if (res_nrm <= res_tol) {
		// mexPrintf("CUDA CG exited early w/ residual norm %e < %e\n", res_nrm, res_tol);
		tmp_vec.free();
		res.free();
		prev_res.free();

		CHECK_CUSPARSE(cusparseDestroySpMat(mat_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(rhs_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(sol_cuda));
		CHECK_CUSPARSE(cusparseDestroyDnVec(tmp_cuda));
		return sol;
	}

	DevPtr<double> p            = new_gpu_buffer<double>(res.size);
	cusparseDnVecDescr_t p_cuda = cuda_vec(p);
	CHECK_CUBLAS(gpu_math.copy(p, res));

	uint32_t current_iter = 0;
	while (current_iter < maxnit) {
		current_iter += 1;
		CHECK_CUSPARSE(gpu_math.sparse_matvec(mat_cuda, p_cuda, tmp_cuda));
		double alpha = res_nrm * res_nrm / CHECK_CUBLAS_VARIANT(gpu_math.dot_product(tmp_vec, p));

		// Update solution. Calculate x_k+1 = x_k + ap
		CHECK_CUBLAS(gpu_math.axpy(alpha, p, sol));

		// Update residual.
		CHECK_CUBLAS(gpu_math.copy(prev_res, res));
		CHECK_CUBLAS(gpu_math.axpy(-alpha, tmp_vec, res));
		double prev_res_nrm = res_nrm;
		res_nrm             = CHECK_CUBLAS_VARIANT(gpu_math.norm(res));

		if (res_nrm < res_tol) { break; }

		double beta = res_nrm * res_nrm / prev_res_nrm / prev_res_nrm;

		CHECK_CUBLAS(gpu_math.scale(p, beta));
		CHECK_CUBLAS(gpu_math.axpy(1, res, p));
	}
#pragma endregion

	// mexPrintf("CUDA CG stopped after %d iterations w/ res norm %e.\n", current_iter, res_nrm);

	tmp_vec.free();
	res.free();
	prev_res.free();

	CHECK_CUSPARSE(cusparseDestroySpMat(mat_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(rhs_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(sol_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(tmp_cuda));
	CHECK_CUSPARSE(cusparseDestroyDnVec(p_cuda));

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
	if (nrhs < 2) { mexErrMsgIdAndTxt("Matlode:cudaCG:nrhs", "Two inputs required"); }

	if (nlhs != 1) { mexErrMsgIdAndTxt("Matlode:cudaCG:nlhs", "One output required"); }

	// First Argument - A (m-by-n matrix of doubles)
	// TODO: make work with floats as well

	const size_t* a_dimensions = mxGetDimensions(prhs[0]);
	size_t a_num_dimension     = mxGetNumberOfDimensions(prhs[0]);
	if (a_num_dimension != 2) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:a_dims", "Input 1 (A) must be n-by-n dimension (matrix)");
	}

	if (a_dimensions[0] != a_dimensions[1]) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:a_dims", "Input 1 (A) must be square matrix");
	}

	if (!mxIsSparse(prhs[0])) { mexErrMsgIdAndTxt("Matlode:cudaCG:a_sparse", "Input 1 (A) must be sparse"); }

	// Second Argument - b (n-row vector of doubles)

	if (mxGetNumberOfDimensions(prhs[1]) > 2 ||
	    mxGetNumberOfDimensions(prhs[1]) == 2 &&
	        (mxGetDimensions(prhs[1])[0] != 1 && mxGetDimensions(prhs[1])[1] != 1)) {
		std::ostringstream ss;
		ss << "Input 2 (b) must be a vector. Instead found ";
		size_t i;
		for (i = 0; i < mxGetNumberOfDimensions(prhs[1]) - 1; i++) { ss << mxGetDimensions(prhs[1])[i] << 'x'; }
		ss << mxGetDimensions(prhs[1])[i];
		mexErrMsgIdAndTxt("Matlode:cudaCG:b_dim", ss.str().c_str());
	}

	if (mxGetNumberOfElements(prhs[1]) != a_dimensions[1]) {
		std::ostringstream ss;
		ss << "Input 2 (b) length must match Input 1 (A) matrix dimension. " << mxGetNumberOfElements(prhs[1])
		   << " v.s. " << a_dimensions[0] << 'x' << a_dimensions[1];
		mexErrMsgIdAndTxt("Matlode:cudaCG:b_dim", ss.str().c_str());
	}

	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:b_type", "Input 2 (b) must be noncomplex double");
	}

	// Optional arguments

	if (nrhs < 3) { return; }

	// Third argument - tolerance

	if (mxGetNumberOfElements(prhs[2]) != 1) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:tol_dim", "Input 3 (tolerance) must be a scalar");
	}

	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:tol_type", "Input 3 (tolerance) must be noncomplex double");
	}

	if (mxGetDoubles(prhs[2])[0] <= 0) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:tol_pos", "Input 3 (tolerance) must be positive");
	}

	if (nrhs < 4) { return; }
	// Fifth argument - maximum iterations

	if (mxGetNumberOfElements(prhs[3]) != 1) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:maxit_dim", "Input 4 (maxit) must be a scalar");
	}

	if (!get_uint(prhs[3]).has_value()) {
		mexErrMsgIdAndTxt("Matlode:cudaCG:maxit_type", "Input 4 (maxit) must be convertible to Uint32");
	}
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	checkArguments(nlhs, plhs, nrhs, prhs);

	cudaSetDevice(0);

	// matrix A
	SparseMatCSC<HostPtr> a_host = get_csc(prhs[0]);
	SparseMatCSC<DevPtr> a_dev   = move_csc(a_host);

	// Vector b
	HostPtr<double> b_host = {mxGetDoubles(prhs[1]), a_host.cols, true};
	DevPtr<double> b_dev   = new_gpu_buffer(b_host);

	double tol = 1e-6;
	if (nrhs >= 3) { tol = mxGetDoubles(prhs[2])[0]; }

	uint32_t maxit = std::max(a_host.cols > UINT32_MAX ? UINT32_MAX : (uint32_t) a_host.cols, 10u);
	if (nrhs >= 4) { maxit = get_uint(prhs[3]).value_or(maxit); }

	DevPtr<double> sol_dev = conj_grad(maxit, tol, a_dev, b_dev);

	// Output vector
	plhs[0]             = mxCreateDoubleMatrix(a_host.rows, 1, mxREAL);
	HostPtr<double> out = {mxGetDoubles(plhs[0]), a_host.rows, true};

	CHECK_CUDA(cuda_memcpy(out, sol_dev));

	b_dev.free();
	sol_dev.free();
}