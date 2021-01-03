#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>

#include "matrix.h"

const size_t N = 40;
const size_t Width = 92;
const size_t Height = 112;
const size_t Eigenfaces = 28;
const std::string Path = "data/";
const size_t Res = Width * Height;
const std::string SampleName = "10";
const size_t MaxValue = 255;

std::pair<std::vector<double>, Matrix> eigenface(Matrix *mat);
