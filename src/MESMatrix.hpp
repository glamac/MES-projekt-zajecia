#ifndef _MES_MATRIX_H
#define _MES_MATRIX_H

#include <cassert>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <ostream>

namespace MES {
    typedef double real_t;
    class Matrix {
        real_t* data;
        size_t rows;
        size_t cols;
    public:
        Matrix(size_t rows = 4, size_t cols = 4, real_t initial_value = 0);
        Matrix(size_t rows, size_t cols, std::initializer_list<real_t> values);
        Matrix(const Matrix& other); // Copy constructor
        Matrix(Matrix&& other) noexcept; // Move constructor
        ~Matrix(); // Destructor

        // Assignment
        Matrix& operator=(const Matrix& other); // Copy assignment
        Matrix& operator=(Matrix&& other) noexcept; // Move assignment
        Matrix& operator=(std::initializer_list<real_t> values); // Assign initializer list
        // Access
        const real_t& operator[](size_t row, size_t col = 0) const;
        real_t& operator[](size_t row, size_t col = 0);

        // Operations
        Matrix transpose() const;
        Matrix operator-() const;
        Matrix& operator+=(const Matrix& rhs);
        Matrix& operator-=(const Matrix& rhs);
        Matrix& operator*=(real_t alpha);
        Matrix& operator/=(real_t alpha);

        // Utility
        real_t min() const;
        real_t max() const;

        size_t numRows() const {return rows;}
        size_t numCols() const {return cols;}

        friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
        friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
        friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
        friend Matrix operator*(const Matrix& lhs, real_t alpha);
        friend Matrix operator*(real_t alpha, const Matrix& rhs);
        friend Matrix operator/(const Matrix& lhs, real_t alpha);
        friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);

    };
    Matrix gauss(const Matrix& A, const Matrix& B);

}
#endif // _MES_MATRIX_H
