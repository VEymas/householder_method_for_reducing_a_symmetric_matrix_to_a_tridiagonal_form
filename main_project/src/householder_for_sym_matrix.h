#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <cstdlib>

namespace {
    const double EPS = std::numeric_limits<double>::epsilon();
}

namespace householder {

//res m*k A m*n B n*k
std::vector<double> matrix_multiplication (const std::vector<double>& matrix1, const std::vector<double>& matrix2, int m, int n, int k);

std::vector<double> mul_matrix_by_number(const std::vector<double>& matrix, double num);

double vector_norm(const std::vector<double>& vector);

double vector_norm2(const std::vector<double>& vector);

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size);

std::vector<double> householder_method(const std::vector<double>& matrix_, int size);
}