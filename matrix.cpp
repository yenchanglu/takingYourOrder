#include "matrix.h"
#include <iomanip>

Matrix::Matrix() {
    nrow = 0;
    ncol = 0;
};

Matrix::Matrix(size_t num_row, size_t num_col) {
    nrow = num_row;
    ncol = num_col;
    for (size_t r = 0; r < nrow; ++r) {
        std::vector<double> row (ncol, 0);
        array.push_back(row);
    }
}

Matrix::Matrix(size_t num_row, size_t num_col, std::vector<std::vector<double> > elements) {
    nrow = num_row;
    ncol = num_col;
    if (elements.size() != num_row) {
        std::cout << "Warning: number of nrow does not match elements." << std::endl;
    }
    if (nrow > 0 && elements[0].size() != num_col) {
        std::cout << "Warning: number of ncol does not match elements." << std::endl;
    }
    array = elements;
}

Matrix Matrix::transpose() {
    std::vector<std::vector<double> > ret;
    for (size_t r = 0; r < ncol; ++r) {
        std::vector<double> row;
        for (size_t c = 0; c < nrow; ++c) {
            row.push_back(array[c][r]);
        }
        ret.push_back(row);
    }

    return Matrix(ncol, nrow, ret);
}

Matrix Matrix::get_row(size_t num_row) const {
    if (num_row < 0 || num_row >= nrow) {
        throw std::out_of_range("Row number is out of range");
    }

    std::vector<std::vector<double> > row;
    row.push_back(array[num_row]);

    return Matrix(1, ncol, row);
}

Matrix operator* (const Matrix &mat1, const Matrix &mat2) {
    if (mat1.ncol != mat2.nrow) {
    	throw std::out_of_range("Dimesions don't match");
    }

    Matrix ret = Matrix(mat1.nrow, mat2.ncol);
    for (size_t r = 0; r < mat1.nrow; r++) {
        for (size_t c = 0; c < mat2.ncol; ++c) {
            for (size_t k = 0; k < mat1.ncol; ++k) {
                ret.array[r][c] += mat1.array[r][k] * mat2.array[k][c];
            }
        }
    }
    return ret;
}