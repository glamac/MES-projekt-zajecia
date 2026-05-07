#include "MESMatrix.hpp"
#include <algorithm>

namespace MES {
    Matrix::Matrix(size_t rows, size_t cols, real_t initial_value) :
        rows(rows), cols(cols)
    {
        data = new real_t[rows * cols];
        std::fill(data, data + rows*cols, initial_value);
    }

    Matrix::Matrix(size_t rows, size_t cols, std::initializer_list<real_t> values) :
        rows(rows), cols(cols)
    {
        assert(values.size() == 1 || values.size() == rows*cols);
        data = new real_t[rows * cols];
        if(values.size() == 1) {
            std::fill(data, data + rows*cols, *values.begin());
        } else {
            std::copy(values.begin(), values.end(), data);
        }
    }

    Matrix::Matrix(const Matrix& other) : // Copy constructor
        rows(other.rows), cols(other.cols)
    {
        data = new real_t[rows * cols];
        std::copy(other.data, other.data + rows * cols, data);
    }

    Matrix::Matrix(Matrix&& other) noexcept : // Move constructor
        rows(other.rows), cols(other.cols), data(other.data)
    {
        other.data = nullptr;
    }

    Matrix::~Matrix() {
        delete[] data;
    }
    Matrix& Matrix::operator=(const Matrix& other) { // Copy assignment
        if (this == &other) goto ret;
        delete[] data;
        rows = other.rows;
        cols = other.cols;
        data = new real_t[rows * cols];
        std::copy(other.data, other.data + rows*cols, data);
    ret:
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& other) noexcept {
        if (this == &other) goto ret;
        delete[] data;
        rows = other.rows;
        cols = other.cols;
        data = other.data;
        other.rows = 0; other.cols = 0;
        other.data = nullptr;
    ret:
        return *this;
    }

    Matrix& Matrix::operator=(std::initializer_list<real_t> values) {
        assert(values.size() == 1 || values.size() == rows*cols);
        if(values.size() == 1) {
            std::fill(data, data + rows*cols, *values.begin());
        } else {
            std::copy(values.begin(), values.end(), data);
        }
        return *this;
    }

    const real_t& Matrix::operator[](size_t row, size_t col) const {
        return data[row * cols + col];
    }

    real_t& Matrix::operator[](size_t row, size_t col) {
        return data[row * cols + col];
    }

    Matrix Matrix::transpose() const {
        Matrix ret(numCols(), numRows());
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                ret[j, i] = (*this)[i, j];
            }
        }
        return ret;
    }

    Matrix Matrix::operator-() const {
        Matrix ret(numRows(), numCols());
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                ret[i, j] = (*this)[i, j] * (-1);
            }
        }
        return ret;
    }

    Matrix& Matrix::operator+=(const Matrix& rhs) {
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                (*this)[i, j] += rhs[i, j];
            }
        }
        return *this;
    }

    Matrix& Matrix::operator-=(const Matrix& rhs) {
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                (*this)[i, j] -= rhs[i, j];
            }
        }
        return *this;
    }

    Matrix& Matrix::operator*=(const real_t alpha) {
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                (*this)[i, j] *= alpha;
            }
        }
        return *this;
    }

    Matrix& Matrix::operator/=(const real_t alpha) {
        for(size_t i = 0; i < numRows(); i++) {
            for(size_t j = 0; j < numCols(); j++) {
                (*this)[i, j] /= alpha;
            }
        }
        return *this;
    }

    // Utility

    real_t Matrix::min() const {
    	real_t min = (*this)[0,0];
    	for(size_t i = 0; i < numRows(); i++) {
     		for(size_t j = 0; j < numCols(); j++) {
       			if((*this)[i, j] < min) min = (*this)[i, j];
       		}
     	}
     	return min;
    }

    real_t Matrix::max() const {
    	real_t max = (*this)[0,0];
    	for(size_t i = 0; i < numRows(); i++) {
     		for(size_t j = 0; j < numCols(); j++) {
       			if((*this)[i, j] > max) max = (*this)[i, j];
       		}
     	}
     	return max;
    }

    Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
	    assert(lhs.numRows() == rhs.numRows());
		assert(lhs.numCols() == rhs.numCols());
	    Matrix ret(lhs.numRows(), lhs.numCols());
	    for(size_t i = 0; i < lhs.numRows(); i++) {
	        for(size_t j = 0; j < rhs.numCols(); j++) {
	            ret[i, j] = lhs[i, j] + rhs[i, j];
	        }
	    }
	    return ret;
    }

    Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
        return (lhs + (-rhs));
    }

    Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
        assert(lhs.numCols() == rhs.numRows());
        Matrix ret(lhs.numRows(), rhs.numCols());
        for(size_t i = 0; i < lhs.numRows(); i++) {
            for(size_t j = 0; j < rhs.numCols(); j++) {
                real_t sum = 0.0;
                for(size_t k = 0; k < lhs.numCols(); k++) {
                    sum += lhs[i, k] * rhs[k, j];
                }
                ret[i, j] = sum;
            }
        }
        return ret;
    }

    Matrix operator*(const Matrix& lhs, real_t alpha) {
        Matrix ret(lhs.numRows(), lhs.numCols());
        for(size_t i = 0; i < lhs.numRows(); i++) {
            for(size_t j = 0; j < lhs.numCols(); j++) {
                ret[i, j] = lhs[i, j] * alpha;
            }
        }
        return ret;
    }

    Matrix operator*(real_t alpha, const Matrix& rhs) {
        return (rhs * alpha);
    }

    Matrix operator/(const Matrix& lhs, real_t alpha) {
        Matrix ret(lhs.numRows(), lhs.numCols());
        for(size_t i = 0; i < lhs.numRows(); i++) {
            for(size_t j = 0; j < lhs.numCols(); j++) {
                ret[i, j] = lhs[i, j] / alpha;
            }
        }
        return ret;
    }

    std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
        os << "matrix(rows: " << mat.numRows() << ", cols: " << mat.numCols() <<"):\n";
        for(size_t i = 0; i < mat.numRows(); i++) {
            for(size_t j = 0; j < mat.numCols(); j++) {
                os << mat[i, j] << ' ';
            }
            os << '\n';
        }
        return os;
    }

    /// Solves a system of equations
    /// A * X = B; solves for X
    Matrix gauss(const Matrix& A, const Matrix& B) {
        assert(A.numCols() == A.numRows() && A.numRows() == B.numRows() && B.numCols() == 1);
        const size_t N = A.numRows();
        Matrix AB(N, N+1);
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                AB[i, j] = A[i, j];
            }
        }
        for(int i = 0; i < N; i++) {
            AB[i, N] = B[i, 0];
        }
        double m = 0, s = 0;
        for(int i = 0; i < N-1; i++) {
            if(std::abs(AB[i, i]) < 1e-12) {
                std::cout << AB << std::endl << "failed at " << i << " (" << AB[i, i] << ")\n";
                return Matrix(N, 1, -2137);
            }
            for(int j = i+1; j < N; j++) {
                m = -AB[j, i] / AB[i, i];
                for(int k = i+1; k <= N; k++) {
                    AB[j, k] += m * AB[i, k];
                }
            }
        }
        Matrix X(N, 1);
        for(int i = N-1; i >= 0; i--) {
            s = AB[i, N];
            for(int j = N-1; j > i; j--) {
                s -= AB[i, j] * X[j, 0];
            }
            X[i, 0] = s / AB[i, i];
        }
        return X;
    }

}
