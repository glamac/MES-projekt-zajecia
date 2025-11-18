#ifndef _MES_MATRIX_H
#define _MES_MATRIX_H
#include <cassert>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <ostream>

#include "MESCore.hpp"

namespace MES {
    template<size_t matRows = 4, size_t matCols = matRows>
    struct Matrix {
        real_t matrix[matRows][matCols];
        constexpr Matrix() {
            for(int i = 0; i < matRows; i++) {
                for(int j = 0; j < matCols; j++) {
                    matrix[i][j] = 0;
                }
            }
        }

        constexpr Matrix(std::initializer_list<MES::real_t> A) {
            if(A.size() == 1) {
                memset(matrix, A.begin()[0], sizeof matrix);
                return;
            }
            for(int i = 0; i < matRows; i++) {
                for(int j = 0; j < matCols; j++) {
                    matrix[i][j] = A.begin()[j + matCols*i];
                }
            }
        }
        
        constexpr decltype(auto) operator[](this auto& self, size_t row, size_t col) {
            assert(row < matRows);
            assert(col < matCols);
            return self.matrix[row][col];
        }

        Matrix<matCols, matRows> transpose() {
            Matrix<matCols, matRows> T;
            for(int i = 0; i < matRows; i++) {
                for(int j = 0; j < matCols; j++) {
                    T[j, i] = matrix[i][j];
                }
            }
            return T;
        }

        Matrix<matRows, matCols>& operator*=(real_t alpha) {
            for(int i = 0; i < matRows; i++) {
                for(int j = 0; j < matCols; j++) {
                    matrix[i][j] *= alpha;
                }
            }
            return *this;
        }

        friend std::ostream& operator<<(std::ostream &stream, Matrix<matRows, matCols> mat) {
            std::cout << "matrix<" << matRows << ", " << matCols << ">:\n";
            for(int i = 0; i < matRows; i++) {
                for(int j = 0; j < matCols; j++) {
                    std::cout << mat[i, j] << ' ';
                }
                std::cout << std::endl;
            }
            return stream;
        }
    };
    // Matrix * Matrix
    template<size_t lhsRows, size_t lhsCols, size_t rhsCols>
    constexpr Matrix<lhsRows, rhsCols> operator*(
          Matrix<lhsRows, lhsCols> lhs,
          const Matrix<lhsCols, rhsCols >& rhs
        ) {
        Matrix<lhsRows, rhsCols> ret;
        for(int i = 0; i < lhsRows; i++) {
          for(int j = 0; j < rhsCols; j++) {
            ret[i, j] = 0;
          }
        }
        for(int i = 0; i < lhsRows; i++) {
            for(int j = 0; j < rhsCols; j++) {
                for(int k = 0; k < lhsCols; k++) {
                    ret[i, j] += lhs[i, k] * rhs[k, j];
                }
            }
        }
        return ret;
    }

    // Matrix * double
    template<size_t rows, size_t cols>
    constexpr Matrix<rows, cols> operator*(
        Matrix<rows, cols> lhs,
        const real_t &alpha
    ) {
        Matrix<rows, cols> ret;
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                ret[i, j] = lhs[i, j] * alpha;
            }
        }
        return ret;
    }
    template<size_t rows, size_t cols>
    constexpr Matrix<rows, cols> operator*(
        real_t alpha,
        const Matrix<rows, cols> &rhs
    ) { return rhs * alpha; }


    template<size_t rows, size_t cols>
    constexpr Matrix<rows, cols> operator+(
        Matrix<rows, cols> lhs,
        const Matrix<rows, cols> &rhs
    ) {
        Matrix<rows, cols> ret;
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                ret[i, j] = lhs[i, j] + rhs[i, j];
            }
        }
        return ret;
    }

    template<size_t rows, size_t cols>
    constexpr bool operator==(
        Matrix<rows, cols> lhs,
        const Matrix<rows, cols> &rhs
    ) {
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                if(lhs[i, j] != rhs[i, j])
                    return false;
            }
        }
        return true;
    }

}

#endif // _MES_MATRIX_H