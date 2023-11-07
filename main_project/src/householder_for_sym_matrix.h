#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

namespace {
    const double EPS = std::numeric_limits<double>::epsilon();
}

std::vector<double> get_right_down_block(const std::vector<double>& matrix, int size, int step);

std::vector<double> get_left_down_block(const std::vector<double>& matrix, int size, int step);

std::vector<double> get_left_top_block(const std::vector<double>& matrix, int size, int step);

std::vector<double> get_right_top_block(const std::vector<double>& matrix, int size, int step);

std::vector<double> get_matrix_from_blocks(const std::vector<double>& left_top_block, const std::vector<double>& right_top_block,
                        const std::vector<double>& left_down_block, const std::vector<double>& right_down_block, int size, int step);

//res m*k A m*n B n*k
std::vector<double> matrix_multiplication (const std::vector<double>& matrix1, const std::vector<double>& matrix2, int m, int n, int k);

std::vector<double> mul_matrix_by_number(const std::vector<double>& matrix, double num);

std::vector<double> matrix_subtract(const std::vector<double>& matrix1, const std::vector<double>& matrix2);

std::vector<double> matrix_subtract_with_init_u(const std::vector<double>& matrix1, const std::vector<double>& matrix2, int size,
                        int step, std::vector<double>& vector_u);

double vector_norm(const std::vector<double>& vector);

double vector_norm2(const std::vector<double>& vector);

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size);

std::vector<double> householder_method(const std::vector<double>& matrix_, int size);