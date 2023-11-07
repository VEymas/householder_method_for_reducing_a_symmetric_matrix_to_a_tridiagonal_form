#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

namespace {
    const double EPS = std::numeric_limits<double>::epsilon();
}

std::vector<double> get_right_down_block(const std::vector<double>&, int, int);

std::vector<double> get_left_down_block(const std::vector<double>&, int, int);

std::vector<double> get_left_top_block(const std::vector<double>&, int, int);

std::vector<double> get_right_top_block(const std::vector<double>&, int, int);

std::vector<double> get_matrix_from_blocks(const std::vector<double>&, const std::vector<double>&,
                        const std::vector<double>&, const std::vector<double>&, int, int);

std::vector<double> matrix_multiplication (const std::vector<double>&, const std::vector<double>&, int, int, int);

std::vector<double> mul_matrix_by_number(const std::vector<double>&, double);

std::vector<double> matrix_subtraction(const std::vector<double>&, const std::vector<double>&);

std::vector<double> matrix_subtraction_with_init_u(const std::vector<double>&, const std::vector<double>&, int, int, std::vector<double>&);

double vector_norm(const std::vector<double>&);

double vector_norm2(const std::vector<double>&);

void find_the_zeros_of_the_matrix(std::vector<double>&, int);

std::vector<double> householder_method(const std::vector<double>&, int);