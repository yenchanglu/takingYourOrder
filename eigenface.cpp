#include "eigenface.h"

namespace py = pybind11;

Matrix identity_matrix(size_t dimension) {
	Matrix I = Matrix(dimension, dimension);
	for (size_t row = 0; row < dimension; ++row) {
		for (size_t col = 0; col < dimension; ++col) {
			if (row == col) {
				I.array[col][row] = 1;
			}
		}
	}
	return I;
}

std::pair<std::vector<double>, Matrix> sort_by_eigenvalues(std::pair<std::vector<double>, Matrix> eigenface) {
	for (size_t i = 0; i < eigenface.first.size() - 1; ++i) {
		size_t idx = i;
		double val = eigenface.first[i];
		for (size_t  j = i + 1; j < eigenface.first.size(); ++j) {
			if (eigenface.first[j] > val) {
				idx = j;
				val = eigenface.first[j]; 
			}
		}

		if (idx != i) {
			std::swap(eigenface.first[i], eigenface.first[idx]);
			for (size_t row = 0; row < eigenface.second.nrow; ++row) {
				std::swap(eigenface.second.array[row][i], eigenface.second.array[row][idx]);
			}
		}
	}

	return eigenface;
}

std::pair<std::vector<double>, Matrix> eigenface(Matrix *mat) {
	Matrix mat1(mat->nrow, mat->ncol);
	mat1.array = mat->array;

	Matrix mat2 = identity_matrix(mat1.nrow);

	size_t max_rotations = 5 * mat1.nrow * mat1.nrow;

	for (size_t i = 0; i < mat1.nrow; ++i) {
		std::vector<double> row;
		for (size_t j = 0; j < mat1.ncol; ++j) {
			row.push_back(mat->array[i][j]);
		}
		mat1.array.push_back(row);
	}

	std::vector<double> eigenvalues;

	for (size_t it = 0; it < max_rotations; ++it) {
		double max = 0, k, l;
		for (size_t row = 0; row < mat1.nrow - 1; ++row) {
			for (size_t col = row + 1; col < mat1.nrow; ++col) {
				if (fabs(mat1.array[row][col]) >= max) {
					max = fabs(mat1.array[row][col]);
					k = row;
					l = col;
				}
			}
		}
		if (max < 1.0e-12) {
			for (size_t i = 0; i < mat1.nrow; ++i) {
				eigenvalues.push_back(mat1.array[i][i]);
			}
			// Normalize mat2
			for (size_t col = 0; col < mat2.ncol; ++col) {
				double len = 0;
				for (size_t row = 0; row < mat2.nrow; ++row) {
					len += mat2.array[row][col] * mat2.array[row][col];
				}
				for (size_t row = 0; row < mat2.nrow; ++row) {
					mat2.array[row][col] = mat2.array[row][col] / len;
				}
			}

			std::pair<std::vector<double>, Matrix> result(eigenvalues, mat2);
			return sort_by_eigenvalues(result);
		}
		// Perform rotation
		double diff = mat1.array[l][l] - mat1.array[k][k];

		double t;
		if (fabs(mat1.array[k][l]) < fabs(diff) * 1.0e-36) {
			t = mat1.array[k][l] / diff;
		} else {
			double phi = diff / (2.0 * mat1.array[k][l]);
			t = 1.0 / (fabs(phi) + sqrt(phi * phi + 1.0));
			if (phi < 0) {
				t = -t;
			}
		}

		double c = 1.0 / sqrt(t*t + 1.0);
		double s = t * c;
		double tau = s / (1.0 + c);
		double tmp = mat1.array[k][l];
		mat1.array[k][l] = 0;
		mat1.array[k][k] = mat1.array[k][k] - t * tmp;
		mat1.array[l][l] = mat1.array[l][l] + t * tmp;

		for (size_t i = 0; i < k; ++i) {
			tmp = mat1.array[i][k];
			mat1.array[i][k] = tmp - s * (mat1.array[i][l] + tau * tmp);
			mat1.array[i][l] = mat1.array[i][l] + s * (tmp - tau * mat1.array[i][l]);
		}

		for (size_t i = k + 1; i < l; ++i) {
			tmp = mat1.array[k][i];
			mat1.array[k][i] = tmp - s * (mat1.array[i][l] + tau * mat1.array[k][i]);
			mat1.array[i][l] = mat1.array[i][l] + s * (tmp - tau * mat1.array[i][l]);
		}

		for (size_t i = l + 1; i < mat1.nrow; ++i) {
			tmp = mat1.array[k][i];
			mat1.array[k][i] = tmp - s * (mat1.array[l][i] + tau * tmp);
			mat1.array[l][i] = mat1.array[l][i] + s * (tmp - tau * mat1.array[l][i]);
		}

		for (size_t i = 0; i < mat1.nrow; ++i) {
			tmp = mat2.array[i][k];
			mat2.array[i][k] = tmp - s * (mat2.array[i][l] + tau * mat2.array[i][k]);
			mat2.array[i][l] = mat2.array[i][l] + s * (tmp - tau * mat2.array[i][l]);
		}
	}

	// Jacobi method didn't converge
	for (size_t i = 0; i < mat1.nrow; ++i) {
		eigenvalues.push_back(mat1.array[i][i]);
	}
	std::pair<std::vector<double>, Matrix> result(eigenvalues, mat2);

	return sort_by_eigenvalues(result);
}

std::vector<double> read_pgm(std::ifstream &file) {
	std::vector<double> values;
	std::string line;
	getline(file, line);
	getline(file, line);
	getline(file, line);
	getline(file, line);

	size_t val;
	while(file >> val) {
		values.push_back(val);
	}
	return values;
}

void write_pgm(std::string file, Matrix *img) {
	std::stringstream fileName;
	fileName << file;
	std::ofstream image_file(fileName.str().c_str());

	image_file << "P2" << std::endl << Width;
	for (size_t i = 0; i < Res; ++i) {
		size_t val = img->array[0][i];
		if (val < 0) {
			val = 0;
		}
		image_file << val << " ";
	}
	image_file.close();
}

std::vector<std::vector<double> > load_data() {
	std::vector<std::vector<double> > array;
	size_t n_faces = 40;
	size_t n_samples = 9;

	for (size_t face = 0; face < n_faces; ++face) {
		std::vector<std::vector<double> > face_array;
		for (size_t sample = 0; sample < n_samples; ++sample) {
			std::stringstream fileName;
			fileName << Path << "p" << face + 1 << "/" << sample + 1 << ".pgm";
			std::ifstream image(fileName.str().c_str());

			if (image.is_open()) {
				face_array.push_back(read_pgm(image));
				image.close();
			}
		}

		std::vector<double> mean;
		for (size_t i = 0; i < Res; ++i) {
			double sum = 0;
			for (size_t j = 0; j < n_samples; ++j) {
				sum += face_array[j][i];
			}
			mean.push_back(sum / n_samples);
		}
		array.push_back(mean);
	}

	return array;
}

Matrix scale(Matrix mat, double min = 0, double max = 255) {
	double m_min = mat.array[0][0];
	double m_max = mat.array[0][0];

	for (size_t r = 0; r < mat.nrow; ++r) {
		for (size_t c = 0; c < mat.ncol; ++c) {
			if (mat.array[r][c] < m_min) {
				m_min = mat.array[r][c];
			}
			if (mat.array[r][c] > m_max) {
				m_max = mat.array[r][c];
			}
		}
	}
	double old_range = m_max - m_min;
	double new_range = max - min;

	Matrix ret(mat.nrow, mat.ncol);
	ret.nrow = mat.nrow;
	ret.ncol = mat.ncol;

	for (size_t r = 0; r < mat.nrow; ++r) {
		std::vector<double> row;
		for (size_t c = 0; c < mat.ncol; ++c) {
			row.push_back((mat.array[r][c] - m_min) * (new_range / old_range) + min);
		}
		ret.array.push_back(row);
	}
	return ret;
}

size_t recognition(Matrix X, Matrix M, Matrix U, Matrix W) {
    /* Subtract the mean image */
    for (size_t c = 0; c < Res; ++c) {
        X.array[0][c] -= M.array[0][c];
        if (X.array[0][c] < 0) {
            X.array[0][c] = 0;
        }
    }

    /* Find weights */
    Matrix Wx = Matrix(Eigenfaces, 1);
    for (size_t r = 0; r < Eigenfaces; ++r) {
        Wx.array[r][0] = (U.get_row(r)*X.transpose()).array[0][0];
    }

    /* Find the closest face from the trainig set */
    double min_distance = 0;
    size_t image_number = 0;
    for (size_t image = 0; image < N; ++image) {
        double distance = 0;
        for (size_t eigenface = 0; eigenface < Eigenfaces; ++eigenface) {
            distance += fabs(W.array[eigenface][image] - Wx.array[eigenface][0]);
        }
        if (distance < min_distance || image == 0) {
            min_distance = distance;
            image_number = image;
        }
    }
    return image_number;
}

double validation() {
	Matrix A = Matrix(N, Res, load_data());
    Matrix M = Matrix(1, Res);

    /* Find the mean image */
    for (size_t c = 0; c < Res; ++c) {
        double sum = 0;
        for (size_t r = 0; r < N; ++r) {
            sum += A.array[r][c];
        }
        M.array[0][c] = sum / N;
    }

    /* Output the mean image */
    write_pgm("output/meanimage.pgm", &M);

    /* Subtract the mean from each image */
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < Res; ++c) {
            A.array[r][c] -= M.array[0][c];
            if (A.array[r][c] < 0) {
                A.array[r][c] = 0;
            }
        }
    }

    /* Output the normalized images */
    for (size_t i = 0; i < N; ++i) {
        Matrix image = A.get_row(i);
        std::ostringstream filename;
        filename << "output/normalized/" << i << ".pgm";
        write_pgm(filename.str(), &image);
    }

    /* Find the covariance matrix */
    Matrix C = A*A.transpose();

    /* Find eigenvectors of the covariance matrix */
    Matrix V = eigenface(&C).second.transpose();

    /* Find eigenfaces */
    Matrix U = Matrix(Eigenfaces, Res);
    for (size_t r = 0; r < Eigenfaces; ++r) {
        Matrix eigenface = V.get_row(r) * A;

        U.array[r] = eigenface.array[0];
        double norm = 0;
        for (size_t i = 0; i < U.ncol; i++) {
            norm += pow(U.array[r][i], 2);
        }
        norm = sqrt(norm);
        for (size_t i = 0; i < U.ncol; i++) {
            U.array[r][i] /= norm;
        }
        /* Output eigenface */
        eigenface = scale(U.get_row(r));
        std::ostringstream filename;
        filename << "output/eigenfaces/" << r << ".pgm";
        write_pgm(filename.str(), &eigenface);
    }

    /* Find weights */
    Matrix W = Matrix(Eigenfaces, N);
    for (size_t r = 0; r < Eigenfaces; ++r) {
        for (size_t c = 0; c < N; ++c) {
            W.array[r][c] = (U.get_row(r) * A.get_row(c).transpose()).array[0][0];
        }
    }

    double accuracy = 0;
    for (size_t i = 1; i <= N; ++i) {
        /* Read image */
        std::stringstream filename;
        filename << Path << "p" << i << "/" << SampleName << ".pgm";
        std::ifstream image(filename.str().c_str());
        std::vector< std::vector<double> > array;

        if (image.is_open()) {
            array.push_back(read_pgm(image));
            image.close();
        } else {
            std::cout << "Image was not opened.";
        }
        Matrix X = Matrix(1, Res, array);
        size_t image_number = recognition(X, M, U, W);

        std::cout << i << ". " << image_number + 1 << std::endl;
        if (i == image_number + 1) {
            accuracy = accuracy + 1;
        }
    }

    //std::cout << accuracy / N << std::endl;
    return (accuracy / N);
}

#if 0
int main(int argc, const char * argv[]) {
    /*
    A contains the images as nrow. A is N * Width * Height, [A]i,j is the jth pixel value of the ith image.
    M contains the mean image. M is 1 *M matrix, [M]0,j is the jth pixel value of the mean image.
    C is the covariance matrix (each pixel is a random variable). C is N * N, computed as AA^T.
    V contains eigenvectors of C as nrow. V is N * N.
    U contains eigenfaces as nrow. U is N * Width * Height.
    W contains the weights of the eigenfaces in the training images. W is N * N, [W]i,j is the weight of the ith eigenface in the jth image.
    X is a 1 * M matrix of the image to recognize (normalized).
    Wx is a Nx1 matrix, [W]i,0 is the weight of the ith eigenface in the X.
    */
    Matrix A = Matrix(N, Res, load_data());
    Matrix M = Matrix(1, Res);

    /* Find the mean image */
    for (size_t c = 0; c < Res; ++c) {
        double sum = 0;
        for (size_t r = 0; r < N; ++r) {
            sum += A.array[r][c];
        }
        M.array[0][c] = sum / N;
    }

    /* Output the mean image */
    write_pgm("output/meanimage.pgm", &M);

    /* Subtract the mean from each image */
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < Res; ++c) {
            A.array[r][c] -= M.array[0][c];
            if (A.array[r][c] < 0) {
                A.array[r][c] = 0;
            }
        }
    }

    /* Output the normalized images */
    for (size_t i = 0; i < N; ++i) {
        Matrix image = A.get_row(i);
        std::ostringstream filename;
        filename << "output/normalized/" << i << ".pgm";
        write_pgm(filename.str(), &image);
    }

    /* Find the covariance matrix */
    Matrix C = A*A.transpose();

    /* Find eigenvectors of the covariance matrix */
    Matrix V = eigenface(&C).second.transpose();

    /* Find eigenfaces */
    Matrix U = Matrix(Eigenfaces, Res);
    for (size_t r = 0; r < Eigenfaces; ++r) {
        Matrix eigenface = V.get_row(r) * A;

        U.array[r] = eigenface.array[0];
        double norm = 0;
        for (size_t i = 0; i < U.ncol; i++) {
            norm += pow(U.array[r][i], 2);
        }
        norm = sqrt(norm);
        for (size_t i = 0; i < U.ncol; i++) {
            U.array[r][i] /= norm;
        }
        /* Output eigenface */
        eigenface = scale(U.get_row(r));
        std::ostringstream filename;
        filename << "output/eigenfaces/" << r << ".pgm";
        write_pgm(filename.str(), &eigenface);
    }

    /* Find weights */
    Matrix W = Matrix(Eigenfaces, N);
    for (size_t r = 0; r < Eigenfaces; ++r) {
        for (size_t c = 0; c < N; ++c) {
            W.array[r][c] = (U.get_row(r) * A.get_row(c).transpose()).array[0][0];
        }
    }

    double accuracy = 0;
    for (size_t i = 1; i <= N; ++i) {
        /* Read image */
        std::stringstream filename;
        filename << Path << "p" << i << "/" << SampleName << ".pgm";
        std::ifstream image(filename.str().c_str());
        std::vector< std::vector<double> > array;

        if (image.is_open()) {
            array.push_back(read_pgm(image));
            image.close();
        } else {
            std::cout << "Image was not opened.";
        }
        Matrix X = Matrix(1, Res, array);
        size_t image_number = recognition(X, M, U, W);

        std::cout << i << ". " << image_number + 1 << std::endl;
        if (i == image_number + 1) {
            accuracy = accuracy + 1;
        }
    }

    std::cout << accuracy / N << std::endl;
    return 0;
}
#endif

PYBIND11_MODULE(_FaceRecog, m) {
    m.doc() = "face recognition module";
    m.def("eigenface", &eigenface);
    m.def("load_data", &load_data);
    m.def("read_pgm", &read_pgm);
    m.def("write_pgm", &write_pgm);
    m.def("recognition", &recognition);
    m.def("validation", &validation);

    py::class_<Matrix>(m, "Matrix")        
        .def(py::init<>());
}
