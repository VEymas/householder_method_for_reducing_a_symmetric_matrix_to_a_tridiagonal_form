#include "householder_for_sym_matrix.h"

std::vector<double> get_right_down_block(const std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((size - step - 1) * (size - step - 1));
    int k{-1}, q{0};
    for (int i = step + 1; i < size; ++i) {
        ++k;
        q = 0;
        for (int j = step + 1; j < size; ++j) {
            block[k * (size - step - 1) + q] = matrix[i * size + j];
            ++q;
        }
    }
    return block;
}

std::vector<double> get_left_down_block(const std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((size - 1 - step) * (1 + step));
    int k{-1}, q{0};
    for (int i = 1 + step; i < size; ++i) {
        ++k;
        q = 0;
        for (int j = 0; j < step + 1; ++j) {
            block[k * (step + 1) + q] = matrix[i * size + j];
            ++q;
        }
    }
    return block;
}


std::vector<double> get_left_top_block(const std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((step + 1) * (step + 1));
    for (int i = 0; i < step + 1; ++i) {
        for (int j = 0; j < step + 1; ++j) {
            block[i * (step + 1) + j] = matrix[i * size + j];
        }
    }
    return block;
}

std::vector<double> get_right_top_block(const std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((1 + step) * (size - 1 - step));
    int k{-1}, q{0};
    for (int i = 0; i < step + 1; ++i) {
        ++k;
        q = 0;
        for (int j = 1 + step; j < size; ++j) {
            block[k * (size - 1 - step) + q] = matrix[i * size + j];
            ++q;
        }
    }
    return block;
}

std::vector<double> get_matrix_from_blocks(const std::vector<double>& left_top_block, const std::vector<double>& right_top_block,
                        const std::vector<double>& left_down_block, const std::vector<double>& right_down_block, int size, int step) {
    std::vector<double> matrix(size * size);
    for (int i = 0; i < step + 1; ++i) {
        for (int j = 0; j < step + 1; ++j) {
            matrix[i * size + j] = left_top_block[i * (step + 1) + j];
        }
    }

    int k{-1}, q{0};
    for (int i = 1 + step; i < size; ++i) {
        ++k;
        q = 0;
        for (int j = 0; j < step + 1; ++j) {
            matrix[i * size + j] = left_down_block[k * (step + 1) + q];
            ++q;
        }
    }

    k = -1;
    q = 0;
   for (int i = 0; i < step + 1; ++i) {
        ++k;
        q = 0;
        for (int j = 1 + step; j < size; ++j) {
            matrix[i * size + j] = right_top_block[k * (size - 1 - step) + q];
            ++q;
        }
    }

    k = -1;
    q = 0;
    for (int i = size - (size - step - 1); i < size; ++i) {
        ++k;
        q = 0;
        for (int j = size - (size - step - 1); j < size; ++j) {
             matrix[i * size + j] = right_down_block[k * (size - step - 1) + q];
            ++q;
        }
    }
    return matrix;
}

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

std::vector<double> matrix_subtraction(const std::vector<double>& matrix1, const std::vector<double>& matrix2) {
    if (matrix1.size() != matrix2.size()) {
        std::cerr << "Bad matrix in matrix_subtraction (sizes don't match)" << std::endl;
    }
    int size = matrix1.size();
    std::vector<double> res_matrix(size);
    for (int i = 0; i < size; ++i) {
        res_matrix[i] = matrix1[i] - matrix2[i];
    }
    return res_matrix;
}

std::vector<double> matrix_subtraction_with_init_u(const std::vector<double>& matrix1, const std::vector<double>& matrix2, int size, int step, std::vector<double>& vector_u) {
    if (matrix1.size() != matrix2.size()) {
        std::cerr << "Bad matrix in matrix_subtraction (sizes don't match)" << std::endl;
    }
    int vector_size = (size - step - 1);
    std::vector<double> res_matrix(matrix1.size());
    int j = 0;
    for (int i = 0; i < matrix1.size(); ++i) {
        res_matrix[i] = matrix1[i] - matrix2[i];
        if (i % vector_size == 0 && i != 0) {
            vector_u[j++] = res_matrix[i];
        }
    }

    vector_u.resize(vector_size - 1);
    vector_u[0] += vector_norm(vector_u);
    
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
    std::vector<double> vector_u(size - 1);
    for (int i = 0; i < size - 2; ++i) {
        // init vector_u
        if (i == 0) {
            for (int j = 0; j < size - i - 1; ++j) {
                vector_u[j] = matrix[(j + i + 1) * size + i];
            }
            vector_u[0] += vector_norm(vector_u);
        }
        
        std::vector<double> left_top_block = get_left_top_block(matrix, size, i);
        
        std::vector<double> left_down_block = get_left_down_block(matrix, size, i);
        std::vector<double> tmp_left_down = matrix_multiplication(vector_u, left_down_block, 1, size - 1 - i, 1 + i);
        tmp_left_down = matrix_multiplication(vector_u, tmp_left_down, size - 1 - i, 1, 1  +i);
        tmp_left_down = mul_matrix_by_number(tmp_left_down, 2/vector_norm2(vector_u));
        left_down_block = matrix_subtraction(left_down_block, tmp_left_down);

        std::vector<double> right_top_block = get_right_top_block(matrix, size, i);
        std::vector<double> tmp_right_top = matrix_multiplication(right_top_block, vector_u, 1 + i, size - 1 - i, 1);
        tmp_right_top = matrix_multiplication(tmp_right_top, vector_u, 1 + i, 1, size - 1 - i);
        tmp_right_top = mul_matrix_by_number(tmp_right_top, 2/vector_norm2(vector_u));
        right_top_block = matrix_subtraction(right_top_block, tmp_right_top); 

        std::vector<double> right_down_block = get_right_down_block(matrix, size, i);
        std::vector<double> tmp_right_down = matrix_multiplication(vector_u, right_down_block, 1, size - 1 - i, size - 1 - i);
        tmp_right_down = matrix_multiplication(vector_u, tmp_right_down, size - 1 - i, 1, size - 1 - i);
        tmp_right_down = mul_matrix_by_number(tmp_right_down, 2/vector_norm2(vector_u));
        right_down_block = matrix_subtraction(right_down_block, tmp_right_down);
        tmp_right_down = matrix_multiplication(right_down_block, vector_u, size - 1 - i, size - 1 - i, 1);
        tmp_right_down = matrix_multiplication(tmp_right_down, vector_u, size - 1 - i, 1, size - 1 - i);
        tmp_right_down = mul_matrix_by_number(tmp_right_down, 2/vector_norm2(vector_u));
        right_down_block = matrix_subtraction_with_init_u(right_down_block, tmp_right_down, size, i, vector_u);

        matrix = get_matrix_from_blocks(left_top_block, right_top_block, left_down_block, right_down_block, size, i);
    }
    find_the_zeros_of_the_matrix(matrix, size);
    return matrix;
}