#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <chrono>

namespace {
    const double EPS = std::numeric_limits<float>::epsilon();
}

namespace householder {

void calculate_reflection_vector(int step, int size, double* reflection_vector, double* matrix, double& norm);

void calculate_left(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm);

void calculate_right(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm);

void householder_method(double* matrix, int matrix_size, double* test_reflection_vectors, double* test_norm);

void calculate_error_for_househ(int matrix_size, double* matrix, double* test_reflection_vectors, double* test_norm, double* test_matrix);

}