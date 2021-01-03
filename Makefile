FLAGS ?=

UNAME_S := $(shell uname -s)

CXX = g++
CXXFLAGS := ${CXXFLAGS} \
	-std=c++17 -O3 -g -m64 \
	-Wall -Wextra -shared -fPIC \
	`python3 -m pybind11 --includes` matrix.cpp eigenface.cpp \
	-o _FaceRecog`python3-config --extension-suffix`

.PHONY: module test clean
module:
	$(CXX) $(CXXFLAGS)

clean:
	rm -rf *.so __pycache__ .pytest_cache performance.txt

test: module
	pytest
