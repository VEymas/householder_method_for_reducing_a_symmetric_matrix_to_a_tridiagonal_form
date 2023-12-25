#include "householder_for_sym_matrix.h"

namespace householder{

#ifdef USE_BLAS

extern "C" {
    double dnrm2_(const int *, const double *, const int *);
    void dgemv_(const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
    void dger_(const int *, const int *, const double *, const double *, const int *, const double *, const int *, double *, const int *);
}

void calculate_reflection_vector(int step, int size, double* reflection_vector, double* matrix, double& norm) {
    int ione = 1;
    int len = size - step - 1;
    int pas = step * size;
    for(int j = step + 1; j < size; j++) {
        reflection_vector[j - step - 1] = matrix[pas + j];
    }
    norm = dnrm2_(&len, reflection_vector, &ione);
    reflection_vector[0] -= norm;
    norm *= norm;
    norm -= matrix[pas + step + 1] * matrix[pas + step + 1] - reflection_vector[0] * reflection_vector[0];
    norm /= 2;
}

void calculate_left(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm) {
    char T = 'T';
    int dif = matrix_size - vector_size;
    int i_one = 1, vector_size_1 = vector_size + 1;
    double sum, alpha = 1;
    double* vector_sum = new double[vector_size + 1];
    double rev_norm = -1 / norm, zero = 0;

    dgemv_(&T, &vector_size, &vector_size_1, &rev_norm, matrix + (dif - 1) * matrix_size + dif, &matrix_size, reflection_vector, &i_one, &zero, vector_sum, &i_one);
    dger_(&vector_size, &vector_size_1, &alpha, reflection_vector, &i_one, vector_sum, &i_one, matrix + (dif - 1) * matrix_size + dif, &matrix_size);
    delete[] vector_sum;
}

void calculate_right(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm) {
    char N = 'N';
    double rev_norm = 1 / norm, zero = 0, minus_one = -1;
    int dif = matrix_size - vector_size, i_one = 1;
    int vector_size_1 = vector_size + 1;
    double* vector_sum = new double[vector_size + 1];

    dgemv_(&N, &vector_size_1, &vector_size, &rev_norm, matrix + dif * matrix_size + dif - 1, &matrix_size, reflection_vector, &i_one, &zero, vector_sum, &i_one);
    dger_(&vector_size_1, &vector_size, &minus_one, vector_sum, &i_one, reflection_vector, &i_one,  matrix + dif * matrix_size + dif - 1, &matrix_size);
    delete[] vector_sum;
}

#else

void calculate_reflection_vector(int step, int size, double* reflection_vector, double* matrix, double& norm) {
    norm = 0;
    int pas = step * size;
    for(int j = step + 1; j < size; j++) {
        reflection_vector[j - step - 1] = matrix[pas + j];
        norm += matrix[pas + j] * matrix[pas + j];
    }
    if (reflection_vector[0] > 0) {
        reflection_vector[0] += std::sqrt(norm);
    } else {
        reflection_vector[0] -= std::sqrt(norm);
    }
    norm -= matrix[pas + step + 1] * matrix[pas + step + 1] - reflection_vector[0] * reflection_vector[0];
    norm /= 2;
}

void calculate_left(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm) {
    int dif = matrix_size - vector_size;
    double sum;
    for(int i = dif - 1; i < matrix_size; i++) {
        sum = 0;
        for(int j = dif; j < matrix_size; j++) {
            sum += (matrix[i *  matrix_size + j] * reflection_vector[j - dif]);
        }
        sum /= norm;
        for(int j = dif; j < matrix_size; j++) {
            matrix[i *  matrix_size + j] -= (reflection_vector[j - dif] * sum);
        }
    }
}

void calculate_right(double* matrix, double* reflection_vector, int matrix_size, int vector_size, double norm) {
int dif = matrix_size - vector_size;
    
    double* vector_sum = new double[matrix_size - dif + 1];

    for (int i = 0; i < matrix_size - dif + 1; ++i) {
        vector_sum[i] = 0;
    }

    for(int j = dif; j < matrix_size; j++) {
        for(int i = dif - 1; i < matrix_size; i++) {
            vector_sum[i - dif + 1] += (matrix[j * matrix_size + i] * reflection_vector[j - dif]);
        }
    }

    for(int i = 0; i < matrix_size - dif + 1; i++) {
        vector_sum[i] /= norm;
    }

    for(int j = dif; j < matrix_size; j++) {
        for(int i = dif - 1; i < matrix_size; i++) {
            matrix[j * matrix_size + i] -= (reflection_vector[j - dif] * vector_sum[i - dif + 1]);
        }
    }
}

#endif

void householder_method(double* matrix, int matrix_size, double* test_reflection_vectors, double* test_norm) {   
    double norm = 0;
    double* x = new double[matrix_size - 1];
    for(int i = 0; i < matrix_size - 2; i++) {
        calculate_reflection_vector(i, matrix_size, x, matrix, norm);
        calculate_left(matrix, x, matrix_size, matrix_size - i - 1, norm);
        calculate_right(matrix, x, matrix_size, matrix_size - i - 1, norm);
        for(int j = 0; j < matrix_size - i - 1; j++) {
            test_reflection_vectors[i * (matrix_size - 1) + j] = x[j];
        }
        test_norm[i] = norm;
    }
}

void calculate_error_for_househ(int matrix_size, double* matrix, double* test_reflection_vectors, double* test_norm, const double* test_matrix) {
    double error = 0, matrix_norm = 0;
    for(int i = matrix_size - 3; i >= 0; i--) {
        calculate_left(matrix, test_reflection_vectors + i * (matrix_size - 1), matrix_size, matrix_size - i - 1, test_norm[i]);
        calculate_right(matrix, test_reflection_vectors + i * (matrix_size - 1), matrix_size, matrix_size - i - 1, test_norm[i]);
    }

    for(int i = 0; i < matrix_size; i++) {
        for(int j = 0; j < matrix_size; j++) {
            matrix_norm += test_matrix[i * matrix_size + j] * test_matrix[i * matrix_size + j];
            error += (test_matrix[i * matrix_size + j] - matrix[i * matrix_size + j]) * (test_matrix[i * matrix_size + j] - matrix[i * matrix_size + j]);
        }
    }
    std::cout << "abs error - " << std::sqrt(error) << " rel error - " << std::sqrt(error) / std::sqrt(matrix_norm) << std::endl;
}
} //namespace householder