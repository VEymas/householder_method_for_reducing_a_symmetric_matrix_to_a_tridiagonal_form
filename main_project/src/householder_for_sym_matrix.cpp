#include "householder_for_sym_matrix.h"

namespace householder{

//res m*k A m*n B n*k
std::vector<double> matrix_multiplication (const std::vector<double>& matrix1, const std::vector<double>& matrix2, int m, int n, int k) {
    std::vector<double> res_matrix(m * k, 0);
    for (int i = 0; i < m; ++i) {
        for (int q = 0; q < n; ++q) {
            for (int j = 0; j < k; ++j) {
                res_matrix[i * k + j] += matrix1[i * n + q] * matrix2[q * k + j];
            }
        }
    }
    return res_matrix;
} 

std::vector<double> mul_matrix_by_number(const std::vector<double>& matrix, double num) {
    int size = matrix.size();
    std::vector<double> res_matrix(size);
    for (int i = 0; i < size; ++i) {
        res_matrix[i] = matrix[i] * num;
    }
    return res_matrix;
}

double vector_norm(const std::vector<double>& vector) {
    int n = vector.size();
    double res{};
    for (int i = 0; i < n; ++i) {
        res += vector[i] * vector[i];
    }
    return sqrt(res);
}

double vector_norm2(const std::vector<double>& vector) {
    int n = vector.size();
    double res{};
    for (int i = 0; i < n; ++i) {
        res += vector[i] * vector[i];
    }
    return res;
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

std::vector<double> householder_method(const std::vector<double>& matrix_, int size) {
    std::vector<double> matrix(size * size);
    for (int i = 0; i < size * size; ++i) {
        matrix[i] = matrix_[i];
    }
    std::vector<double> vector_s(size - 1);
    std::vector<double> vector_u;
    for (int i = 0; i < size - 2; ++i) {
        if (i == 0) {
            for (int j = 0; j < size - i - 1; ++j) {
                vector_s[j] = matrix[i * size + j + i + 1];
            }
            vector_u = vector_s;
            vector_u[0] += vector_norm(vector_s);
        }
        for (int j = i + 1; j < size; ++j) {
            if (j == i + 1) {
                matrix[i * size + j] = - vector_norm(vector_s);
                matrix[j * size + i] = - vector_norm(vector_s);
            } else {
                matrix[i * size + j] = 0;
                matrix[j * size + i] = 0;
            }
        }
        double gamma = 2 / vector_norm2(vector_u);
        std::vector<double> p(size - i - 1, 0);
        for (int q = i + 1; q < size; ++q) {
            for (int j = i + 1; j < size; ++j) {
                p[q - i - 1] += matrix[q * size + j] * vector_u[j - i - 1];
            }
        }

        p = mul_matrix_by_number(p, gamma);

        std::vector<double> up = matrix_multiplication(vector_u, p, size - i - 1, 1, size - i - 1);
        std::vector<double> pu = matrix_multiplication(p, vector_u, size - i - 1, 1, size - i - 1);
        std::vector<double> uu = matrix_multiplication(vector_u, vector_u, size - i - 1, 1, size - i - 1);
        double coeff{};
        for (int j = 0; j < vector_u.size(); ++j) {
            coeff += vector_u[j] * p[j];
        }
        coeff *= gamma;
        uu = mul_matrix_by_number(uu, coeff);

        for (int j = 0; j < up.size(); ++j) {
            up[j] += pu[j] - uu[j];
        }

        vector_s.resize(vector_s.size() - 1);
        for (int q = i + 1; q < size; ++q) {
            for (int j = i + 1; j < size; ++j) {
                matrix[q * size + j] -= up[(q - i - 1) * (size - i - 1) + (j - i - 1)];
                if (j == i + 1 && q != i + 1) {
                    vector_s[q - i - 2] = matrix[q * size + j];
                }
            }
        }
        vector_u = vector_s;
        vector_u[0] += vector_norm(vector_s);
    }
    find_the_zeros_of_the_matrix(matrix, size);
    return matrix;
}
}