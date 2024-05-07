#include <optional>
#include <sstream>
#include <variant>

#include "cublas.h"
#include "cusparse.h"
#include "mex.h"

// TODO: Add stack traces
#define CHECK_CUDA(func)                                                                                         \
	{                                                                                                            \
		cudaError_t status = (func);                                                                             \
		if (status != cudaSuccess) {                                                                             \
			std::ostringstream ss;                                                                               \
			ss << "CUDA API failed at line " << __LINE__ << " with error: " << cudaGetErrorString(status) << "(" \
			   << status << ")";                                                                                 \
			mexErrMsgIdAndTxt("Matlode:cudaMex", ss.str().c_str());                                              \
		}                                                                                                        \
	}

#define CHECK_CUSPARSE(func)                                                                                      \
	{                                                                                                             \
		cusparseStatus_t status = (func);                                                                         \
		if (status != CUSPARSE_STATUS_SUCCESS) {                                                                  \
			std::ostringstream ss;                                                                                \
			ss << "CUSPARSE API failed at line " << __LINE__ << " with error: " << cusparseGetErrorString(status) \
			   << "(" << status << ")";                                                                           \
			mexErrMsgIdAndTxt("Matlode:cudaMex", ss.str().c_str());                                               \
		}                                                                                                         \
	}

#define CHECK_CUBLAS(func)                                                                                            \
	{                                                                                                                 \
		cublasStatus_t status = (func);                                                                               \
		if (status != CUBLAS_STATUS_SUCCESS) {                                                                        \
			std::ostringstream ss;                                                                                    \
			ss << "CUBLAS API failed at line " << __LINE__ << " with error: " << cublasGetStatusString(status) << "(" \
			   << status << ")";                                                                                      \
			mexErrMsgIdAndTxt("Matlode:cudaMex", ss.str().c_str());                                                   \
		}                                                                                                             \
	}

#define CHECK_CUBLAS_VARIANT(func)                                                                               \
	[](std::variant<double, cublasStatus_t> x) {                                                                 \
		try {                                                                                                    \
			return std::get<double>(x);                                                                          \
		} catch (const std::bad_variant_access& ex) {                                                            \
			cublasStatus_t status = std::get<cublasStatus_t>(x);                                                 \
			std::ostringstream ss;                                                                               \
			ss << "CUSPARSE API failed at line " << __LINE__ << " with error: " << cublasGetStatusString(status) \
			   << "(" << status << ")";                                                                          \
			mexErrMsgIdAndTxt("Matlode:cudaMex", ss.str().c_str());                                              \
		}                                                                                                        \
	}(func)

// TODO: replace this with "DevBox" which works like an owned pointer. Also add separate view class (or use C++20 span).
/// @brief Wrapper around a device pointr (GPU)
/// @tparam T What type it is pointing to
/// @tparam isView Whether or not the pointer is a view into another buffer (and shouldn't be freed)
template <typename T, bool isView = false>
struct DevPtr {
	/// @brief The pointer
	T* ptr = nullptr;
	/// @brief The size of the buffer
	size_t size = 0;

	DevPtr() = default;

	DevPtr(T* ptr, size_t size) : ptr(ptr), size(size) {}

	DevPtr(const DevPtr<T, isView>& other) = default;
	DevPtr(DevPtr<T, isView>&& other) : ptr(other.ptr), size(other.size) {
		other.ptr  = nullptr;
		other.size = 0;
	}

	DevPtr& operator=(DevPtr&& other) {
		ptr  = other.ptr;
		size = other.size;

		other.ptr  = nullptr;
		other.size = 0;

		return *this;
	}

	/// @brief Free the pointer, using `cudaFree`
	void free() {
		if (ptr != nullptr && !isView) {
			cudaFree(ptr);
			ptr = nullptr;
		}
	}

	/// @brief Extract a view of the current buffer
	/// @param offset The offset from which the view will start
	/// @param sizeThe size of the view
	/// @return A pointer to the beginning of the view
	DevPtr<T, true> get_view(size_t offset, size_t size) {
		mxAssert(offset + size <= DevPtr::size, "View cannot overrun the original buffer");
		return {ptr + offset, size};
	}
};

// TODO: replace this with "HostManagedPtr" which only works with managed memory
/// @brief Wrapper around a host pointer (CPU)
/// @tparam T What type it is pointing to
template <typename T>
struct HostPtr {
	/// @brief The pointer
	T* ptr = nullptr;
	/// @brief The size of the buffer
	size_t size = 0;

	/// @brief Managed pointers are managed by some other memory system (i.e. Matlab), so they shouldn't be freed
	bool managed = true;

	HostPtr() = default;

	HostPtr(T* ptr, size_t size, bool managed) : ptr(ptr), size(size), managed(managed) {}

	HostPtr(const HostPtr<T>& other) = default;
	HostPtr(HostPtr<T>&& other) : ptr(other.ptr), size(other.size), managed(other.managed) {
		other.ptr  = nullptr;
		other.size = 0;
	}

	HostPtr& operator=(HostPtr&& other) {
		ptr     = other.ptr;
		size    = other.size;
		managed = other.managed;

		other.ptr  = nullptr;
		other.size = 0;

		return *this;
	}

	/// @brief Free the pointer, using `delete[]`
	void free() {
		if (ptr != nullptr && !managed) {
			delete[] ptr;
			ptr = nullptr;
		}
	}
};

/// @brief Sparse matrix in CSC - generic over which device it's on
/// @tparam Ptr The backing storage of the buffers - i.e. which device to store the matrix on.
template <template <typename> typename Ptr>
struct SparseMatCSC {
	/// @brief Pointer to an array of the non-zero values of the matrix, stored in column-major order. Length of
	/// exactly `num_non_zero`.
	Ptr<double> vals;

	/// @brief Pointer to an array of the indices of the beginning of each column in `vals`. Length of exactly `cols
	/// + 1`. The last element must be exactly `num_non_zero`, signifying the end of `vals`.
	Ptr<size_t> col_indices;

	/// @brief Pointer to an array of the row indices of the non-zero values of the matrix. Each index corresponds to
	/// the value stored in the same index of `vals`. Length of exactly `num_non_zero`.
	Ptr<size_t> row_indices;

	/// @brief The number of rows, columns in the matrix
	size_t rows, cols;

	/// @brief The number of non-zero elements in the matrix
	size_t num_non_zero;

	SparseMatCSC()                     = default;
	SparseMatCSC(SparseMatCSC&& other) = default;

	~SparseMatCSC() {
		vals.free();
		col_indices.free();
		row_indices.free();
	}
};

/// @brief Copy a buffer on the host into a new buffer on the GPU
/// @tparam T The type that the buffer is storing
/// @param host_buf The host buffer to be copied
/// @return The new GPU buffer
template <typename T>
DevPtr<T> new_gpu_buffer(HostPtr<T> host_buf) {
	DevPtr<T> re = {nullptr, host_buf.size};

	CHECK_CUDA(cudaMalloc(&re.ptr, re.size * sizeof(T)));
	CHECK_CUDA(cudaMemcpy(re.ptr, host_buf.ptr, re.size * sizeof(T), cudaMemcpyHostToDevice));

	return re;
}

/// @brief Create a buffer on the GPU
/// @tparam T The type that the buffer is storing
/// @param size The number of elements (`T`) that the buffer can store
/// @return The new GPU buffer
template <typename T>
DevPtr<T> new_gpu_buffer(size_t size) {
	DevPtr<T> re = {nullptr, size};

	CHECK_CUDA(cudaMalloc(&re.ptr, re.size * sizeof(T)));

	return re;
}

/// @brief Copy data from a device buffer to a host buffer
/// @tparam T The type of data being copied
/// @param[out] dest The host buffer to be copied into
/// @param[in] src The device buffer copied from
/// @return The CUDA error status of the copy
template <typename T>
cudaError_t cuda_memcpy(HostPtr<T> dest, DevPtr<T> src) {
	mxAssert(dest.size >= src.size, "Cannot copy into smaller buffer");
	return cudaMemcpy(dest.ptr, src.ptr, src.size * sizeof(T), cudaMemcpyDeviceToHost);
}

/// @brief Get a CSC Sparse Matrix from a Matlab argument
/// @param in The Matlab argument that is a sparse CSC matrix
/// @return The CSC matrix
SparseMatCSC<HostPtr> get_csc(const mxArray* in) {
	SparseMatCSC<HostPtr> re;

	re.rows = mxGetM(in);
	re.cols = mxGetN(in);

	// Matlab actually lies here - it has potentially allocated extra room and only tells us the "max" number of nonzero
	// elements. So we need to figure the *actual* number of nonzero elements manually. To do that, we need to first
	// load in the col_indices array.
	re.col_indices  = {mxGetJc(in), re.cols + 1, true};
	re.num_non_zero = 0;
	for (size_t i = 0; i < re.cols; i++) { re.num_non_zero += re.col_indices.ptr[i + 1] - re.col_indices.ptr[i]; }

	// Make sure to indicate that these are managed pointers (coming from Matlab)
	// That way they aren't accidentally freed
	re.vals        = {mxGetPr(in), re.num_non_zero, true};
	re.row_indices = {mxGetIr(in), re.num_non_zero, true};

	return re;
}

/// @brief Move a CSC Sparse matrix to the GPU from the CPU
/// @param host The matrix to be mvoed
/// @return The GPU matrix
SparseMatCSC<DevPtr> move_csc(const SparseMatCSC<HostPtr>& host) {
	SparseMatCSC<DevPtr> dev;

	dev.rows         = host.rows;
	dev.cols         = host.cols;
	dev.num_non_zero = host.num_non_zero;

	// Copy data buffers over
	dev.vals        = new_gpu_buffer(host.vals);
	dev.col_indices = new_gpu_buffer(host.col_indices);
	dev.row_indices = new_gpu_buffer(host.row_indices);

	return dev;
}

/// @brief Register a GPU CSC Sparse Matrix with CUDA, allowing it to be used for matrix operations.
/// @param mat The matrix to be registered
/// @return A handle to the registered matrix, which must be destroyed later.
cusparseSpMatDescr_t cuda_csc(const SparseMatCSC<DevPtr>& mat) {
	cusparseSpMatDescr_t re;

	CHECK_CUSPARSE(cusparseCreateCsc(&re, mat.rows, mat.cols, mat.num_non_zero, mat.col_indices.ptr,
	                                 mat.row_indices.ptr, mat.vals.ptr, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I,
	                                 CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

	return re;
}

/// @brief Register a GPU vector with CUDA, allowing it to be used for matrix operations
/// @tparam isView Whather or not the vector is a view of a larger vector
/// @param vec The vector to be registered
/// @return A handle to the registered vector, which must be destroyed later.
template <bool isView>
cusparseDnVecDescr_t cuda_vec(DevPtr<double, isView> vec) {
	cusparseDnVecDescr_t re;
	cusparseCreateDnVec(&re, vec.size, vec.ptr, CUDA_R_64F);
	return re;
}

/// @brief Helper methods for doing sparse math
class CudaSparseMath {
private:
	cusparseHandle_t sparse_handle;
	cublasHandle_t blas_handle;

	// Buffer for extra sparse mat-vec product scratch space
	void* mv_buffer = nullptr;
	size_t mv_buffer_size;

	const double alpha = 1.0;
	const double beta  = 0.0;

public:
	CudaSparseMath() {
		cusparseCreate(&sparse_handle);
		cublasCreate_v2(&blas_handle);
	}

	~CudaSparseMath() {
		cusparseDestroy(sparse_handle);
		cublasDestroy_v2(blas_handle);

		if (mv_buffer != nullptr) cudaFree(mv_buffer);
	}

	/// @brief Prepare a sparse matrix for multiplication. Must be called at least once before `sparse_matvec` for each
	/// matrix.
	/// @param mat A handle to a sparse matrix - see `cuda_csc`
	/// @param x A handle to a vector - see `cuda_vec`
	/// @param y A handle to a vector - see `cuda_vec`
	/// @return CUDA Status
	cusparseStatus_t prepare_sparse_matvec(cusparseSpMatDescr_t mat, cusparseDnVecDescr_t x, cusparseDnVecDescr_t y) {
		if (mv_buffer != nullptr) cudaFree(mv_buffer);

		cusparseStatus_t status =
		    cusparseSpMV_bufferSize(sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat, x, &beta, y,
		                            CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &mv_buffer_size);

		if (status != CUSPARSE_STATUS_SUCCESS) { return status; }

		cudaMalloc(&mv_buffer, mv_buffer_size);

		return cusparseSpMV_preprocess(sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat, x, &beta, y,
		                               CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, mv_buffer);
	}

	/// @brief Matrix multiplication with a sparse matrix. Must call `prepare_sparse_matvec` first. Calculates `y = A*x`
	/// @param mat The matrix to multiply by
	/// @param x[in] Vector to be multiplied
	/// @param y[out] Output vector
	/// @return CUDA Status
	cusparseStatus_t sparse_matvec(cusparseSpMatDescr_t mat, cusparseDnVecDescr_t x, cusparseDnVecDescr_t y) {
		mxAssert(mv_buffer != nullptr, "Must have run prepare_sparse_matvec to prepare for multiplication");

		return cusparseSpMV(sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat, x, &beta, y, CUDA_R_64F,
		                    CUSPARSE_SPMV_ALG_DEFAULT, mv_buffer);
	}

	/// @brief Calculate 2-norm of a vector
	/// @param vec Vector to calculate the 2-norm of
	/// @return A variant containing either a CUDA status (if not succesful) or the calculated norm
	template <bool isView>
	std::variant<double, cublasStatus_t> norm(DevPtr<double, isView> vec) {
		double norm;
		cublasStatus_t status = cublasDnrm2_v2(blas_handle, vec.size, vec.ptr, 1, &norm);

		if (status != CUBLAS_STATUS_SUCCESS) {
			return status;
		} else {
			return norm;
		}
	}

	/// @brief Compute y = a*x + y
	/// @param alpha The scalar multiple
	/// @param x[in] Vector to add
	/// @param y[in,out] Vector to add to
	/// @return CUDA Status
	template <bool isView1, bool isView2>
	cublasStatus_t axpy(double alpha, DevPtr<double, isView1> x, DevPtr<double, isView2> y) {
		mxAssert(x.size == y.size, "Cannot add two differently sized vectors");

		return cublasDaxpy_v2(blas_handle, x.size, &alpha, x.ptr, 1, y.ptr, 1);
	}

	/// @brief Copy `src` to `dest`
	/// @return CUDA Status
	template <bool isView1, bool isView2>
	cublasStatus_t copy(DevPtr<double, isView1> dest, DevPtr<double, isView2> src) {
		mxAssert(dest.size >= src.size, "Cannot copy into smaller buffer");
		return cublasDcopy_v2(blas_handle, src.size, src.ptr, 1, dest.ptr, 1);
	}

	/// @brief Scale a vector i.e. computes x = a*x
	/// @param vec[in,out] Vector to scale
	/// @param scale Scale (scalar multiple)
	/// @return CUDA Status
	template <bool isView>
	cublasStatus_t scale(DevPtr<double, isView> vec, double scale) {
		return cublasDscal_v2(blas_handle, vec.size, &scale, vec.ptr, 1);
	}

	/// @brief Compute dot product of two vectors
	/// @return Dot product if succesful, CUDA status if not
	template <bool isView>
	std::variant<double, cublasStatus_t> dot_product(DevPtr<double, isView> x, DevPtr<double, isView> y) {
		mxAssert(x.size == y.size, "Cannot take dot product of differently sized vectors");
		double prod;
		cublasStatus_t status = cublasDdot_v2(blas_handle, x.size, x.ptr, 1, y.ptr, 1, &prod);

		if (status != CUBLAS_STATUS_SUCCESS) {
			return status;
		} else {
			return prod;
		}
	}
};
