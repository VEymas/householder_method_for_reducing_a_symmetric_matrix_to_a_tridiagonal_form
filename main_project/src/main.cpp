#include "householder_for_sym_matrix.h"
#include <iostream>
#include <ctime>
#include <vector>
#include <string>

void print_matrix(const std::vector<double> matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << matrix[i * n + j] << "\t";
        }
        std::cout << std::endl;
    }
}

std::vector<double> generate_sym_matrix(int n) {
    std::vector<double> matrix(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            matrix[i * n + j] = (double)std::rand() / 10000000;
            matrix[j * n + i] = matrix[i * n + j];
        }
    }
    return matrix;
}

bool is_matrix_tridiagonal(std::vector<double> matrix, int size) {
    for (int i = 2; i < size; ++i) {
        for (int j = 0; j < i - 1; ++j) {
            if (matrix[i * size + j] != 0 || matrix[j * size + i] != 0) {
                return false;
            } 
        }
    }
    return true;
}

int main() {
    //TESTING
    std::srand(std::time(nullptr));

    int n = 512;
    std::vector<double> matrix = generate_sym_matrix(n);
    clock_t start = clock();
    matrix = householder_method(matrix, n);

    clock_t end = clock();
    std::string message;
    if (is_matrix_tridiagonal(matrix, n)) {
        message = "and matrix is tridiagonal :)";
    } else {
        message = "but matrix is not tridiagonal :(";
    }
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << "n = " << n << " TIME = " << seconds << "s "<< message << std::endl;

}