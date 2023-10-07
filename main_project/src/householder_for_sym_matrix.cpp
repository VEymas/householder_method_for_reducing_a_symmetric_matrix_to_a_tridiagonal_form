#include "householder_for_sym_matrix.h"

std::vector<double> matrix_multiplication (std::vector<double>& matrix1, std::vector<double>& matrix2, int n) {
    std::vector<double> res_matrix(n * n);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                res_matrix[i * n + j] += matrix1[i * n + k] * matrix2[k * n + j];
            }
        }
    }
    return res_matrix;
} 

double vector_norm(std::vector<double> vector) {
    int n = vector.size();
    double res{};
    for (int i = 0; i < n; ++i) {
        res += vector[i] * vector[i];
    }
    return sqrt(res);
}

void multiple_matrix_by_a_number(std::vector<double>& matrix, int size, double number) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[i * size + j] *= number;
        }        
    }
}

std::vector<double> generate_unit_matrix(int size) {
    std::vector<double> unit_matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                unit_matrix[i * size + j] = 1;
            } else {
                unit_matrix[i * size + j] = 0;
            }
        }
    }
    return unit_matrix;
}

std::vector<double> generate_matrix_from_u(std::vector<double> u) {
    int size = u.size();
    std::vector<double> matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[i * size + j] = u[i] * u[j];
        }
    }
    return matrix;
}

std::vector<double> matrix_subtraction(std::vector<double> matrix1, std::vector<double> matrix2, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix1[i * size + j] -= matrix2[i * size + j];
        }
    }
    return matrix1;
}

std::vector<double> generate_matrix_u_from_hausholder(std::vector<double> hausholder_matrix, int step, int size) {
    std::vector<double> res_matrix(size * size, 0.0);
    for (int i = 0; i < step + 1; ++i) {
        res_matrix[i * size + i] = 1;
    }
    for (int i = step + 1; i < size; ++i) {
        for (int j = step + 1; j < size; ++j) {
            res_matrix[i * size + j] = hausholder_matrix[(i - step - 1) * (size - step - 1)+ j - step - 1];
        }
    }
    return res_matrix;
}

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (abs(matrix[i * size + j]) < EPS) {
                matrix[i * size + j] = 0;
            }
        }
    }
}

std::vector<double> householder_method(std::vector<double> matrix, int size) {
    std::vector<double> current_housh((size - 1) * (size - 1));
    for (int i = 0; i < size - 2; ++i) {
        std::vector<double> vector_u(size - 1 - i);
        //init vector_u
        for (int j = 0; j < size - i - 1; ++j) {
            vector_u[j] = matrix[(j + i + 1) * size + i];
        }
        //u[0] or u[i], check this later
        vector_u[0] += vector_norm(vector_u);

        std::vector<double> matrix_from_u = generate_matrix_from_u(vector_u);
        multiple_matrix_by_a_number(matrix_from_u, vector_u.size(), 2);
        multiple_matrix_by_a_number(matrix_from_u, vector_u.size(), 1/(vector_norm(vector_u) * vector_norm(vector_u)));
        std::vector<double> hausholder_matrix = matrix_subtraction(generate_unit_matrix(vector_u.size()), matrix_from_u, vector_u.size());

        std::vector<double> matrix_u = generate_matrix_u_from_hausholder(hausholder_matrix, i, size);
        
        matrix = matrix_multiplication(matrix_u, matrix, size);
        matrix = matrix_multiplication(matrix, matrix_u, size);
    }
    find_the_zeros_of_the_matrix(matrix, size);
    return matrix;
}