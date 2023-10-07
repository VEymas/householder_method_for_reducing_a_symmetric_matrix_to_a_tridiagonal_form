#include <vector>
#include <algorithm>
#include <cmath>

namespace {
    const double EPS = 0.0000001;
}

std::vector<double> matrix_multiplication (std::vector<double>& matrix1, std::vector<double>& matrix2, int n);

double vector_norm(std::vector<double> vector);

void multiple_matrix_by_a_number(std::vector<double>& matrix, int size, double number);

std::vector<double> generate_unit_matrix(int size);

std::vector<double> generate_matrix_from_u(std::vector<double> u);

std::vector<double> matrix_subtraction(std::vector<double> matrix1, std::vector<double> matrix2, int size);

std::vector<double> generate_matrix_u_from_hausholder(std::vector<double> hausholder_matrix, int step, int size);

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size);

std::vector<double> householder_method(std::vector<double> matrix, int size);
