#include "householder_for_sym_matrix.h"
#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <cmath>

void print_matrix(const std::vector<double> matrix, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << matrix[i * size + j] << "\t";
        }
        std::cout << std::endl;
    }
}

double* generate_sym_matrix(int n) {
    double* matrix = new double[n * n];
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            matrix[i * n + j] = (double)std::rand() / 10000000;
            matrix[j * n + i] = matrix[i * n + j];
        }
    }
    return matrix;
}

bool is_matrix_tridiagonal(double* matrix, int size) {
    for (int i = 2; i < size; ++i) {
        for (int j = 0; j < i - 1; ++j) {
            if (abs(matrix[i * size + j]) > EPS || abs(matrix[j * size + i]) > EPS) {
                return false;
            } 
        }
    }
    return true;
}

int
main(void)
{
    std::srand(std::time(nullptr));
    int size = 256;
    for (int i = 0; i < 5; ++i) {
        double* matrix = generate_sym_matrix(size);
        double* test_matrix = new double[size * size];
        double* test_reflection_vectors = new double[(size - 1) * (size - 2)];
        double* test_norm = new double[size - 2]; 
        for(int q = 0; q < size; q++) {
            for(int j = 0; j < size; j++) {
                test_matrix[q * size + j] = matrix[q * size + j];
            }
        }
        clock_t start = clock();
        householder::householder_method(matrix, size, test_reflection_vectors, test_norm);
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        std::string message;
        if (is_matrix_tridiagonal(matrix, size)) {
            message = "and matrix is tridiagonal :)";
        } else {
            message = "but matrix is not tridiagonal :(";
        }
        std::cout << "size = " << size << " TIME = " << seconds << "s "<< message << std::endl;
        householder::calculate_error_for_househ(size, matrix, test_reflection_vectors, test_norm, test_matrix);
        size *= 2;
    }
    return 0;
}