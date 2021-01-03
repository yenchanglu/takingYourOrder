#include <iostream>
#include <vector>

class Matrix {
public:
    size_t nrow;
    size_t ncol;
    std::vector<std::vector<double> > array;

    Matrix();
    Matrix(size_t num_row, size_t num_col);
    Matrix(size_t num_row, size_t num_col, std::vector<std::vector<double> > elements);

    Matrix transpose();
    Matrix get_row(size_t number_of_row) const;
    friend Matrix operator* (const Matrix &mat1, const Matrix &mat2);
};
