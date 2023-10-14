#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

namespace {
    const double EPS = std::numeric_limits<double>::epsilon();
}

std::vector<double> get_right_down_block(std::vector<double>& matrix, int size, int needed_size);

std::vector<double> get_left_down_block(std::vector<double>& matrix, int size, int step);

std::vector<double> matrix_multiplication (const std::vector<double>& matrix1, const std::vector<double>& matrix2, int m, int n, int k);

std::vector<double> mul_matrix_by_number(const std::vector<double>& matrix, double num, int size);

std::vector<double> matrix_subtraction(const std::vector<double>& matrix1, const std::vector<double>& matrix2, int size);

double vector_norm(const std::vector<double>& vector);
double vector_norm2(const std::vector<double>& vector);

std::vector<double> generate_household_matrix_from_vecctor_u(std::vector<double>& u);

std::vector<double> generate_matrix_u_from_household(std::vector<double>& hausholder_matrix, int step, int size);

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size);

template <typename T> int sgn(T val);

template <class T> void print(const T& c); 

std::vector<double> householder_method(const std::vector<double>& matrix, int size);

std::vector<double> get_matrix_from_blocks(std::vector<double>& left_top_block, std::vector<double>& right_top_block,
                        std::vector<double>& left_down_block, std::vector<double>& right_down_block, int size, int step);
