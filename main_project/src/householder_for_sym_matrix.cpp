#include "householder_for_sym_matrix.h"

std::vector<double> get_right_down_block(std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((size - step - 1) * (size - step - 1));
    int k{-1}, q{0};
    for (int i = size - (size - step - 1); i < size; ++i) {
        ++k;
        q = 0;
        for (int j = size - (size - step - 1); j < size; ++j) {
            block[k * (size - step - 1) + q] = matrix[i * size + j];
            ++q;
        }
    }
    return block;
}

std::vector<double> get_left_down_block(std::vector<double>& matrix, int size, int step) {
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


std::vector<double> get_left_top_block(std::vector<double>& matrix, int size, int step) {
    std::vector<double> block((step + 1) * (step + 1));
    for (int i = 0; i < step + 1; ++i) {
        for (int j = 0; j < step + 1; ++j) {
            block[i * (step + 1) + j] = matrix[i * size + j];
        }
    }
    return block;
}

std::vector<double> get_matrix_from_blocks(std::vector<double>& left_top_block, std::vector<double>& right_top_block,
                        std::vector<double>& left_down_block, std::vector<double>& right_down_block, int size, int step) {
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
        ++q;
        k = 0;
        for (int j = 1 + step; j < size; ++j) {
            matrix[i * size + j] = right_top_block[k * (step + 1) + q];
            ++k;
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
                res_matrix[i * n + j] += matrix1[i * n + q] * matrix2[q * n + j];
            }
        }
    }
    return res_matrix;
} 

std::vector<double> mul_matrix_by_number(const std::vector<double>& matrix, double num, int size) {
    std::vector<double> res_matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j)
            res_matrix[i * size + j] = matrix[i * size + j] * num;
    }
    return res_matrix;
}

std::vector<double> matrix_subtraction(const std::vector<double>& matrix1, const std::vector<double>& matrix2, int size) {
    std::vector<double> res_matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res_matrix[i * size + j] = matrix1[i * size + j] - matrix2[i * size + j];
        }
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

std::vector<double> generate_household_matrix_from_vector_u(std::vector<double>& u) {
    int size = u.size();
    double norm = vector_norm(u);
    std::vector<double> matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                matrix[i * size + j] = 1 - u[i] * u[j] * 2 / (norm * norm);
            } else {
                matrix[i * size + j] = -u[i] * u[j] * 2 / (norm * norm);
            }

        }
    }
    return matrix;
}

std::vector<double> generate_matrix_u_from_household(std::vector<double>& householder_matrix, int step, int size) {
    std::vector<double> res_matrix(size * size, 0.0);
    for (int i = 0; i < step + 1; ++i) {
        res_matrix[i * size + i] = 1;
    }
    for (int i = step + 1; i < size; ++i) {
        for (int j = step + 1; j < size; ++j) {
            res_matrix[i * size + j] = householder_matrix[(i - step - 1) * (size - step - 1)+ j - step - 1];
        }
    }
    return res_matrix;
}

void find_the_zeros_of_the_matrix(std::vector<double>& matrix, int size) {
    //std::cout << size << ' ' << matrix.size() << std::endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (abs(matrix[i * size + j]) < EPS) {
                matrix[i * size + j] = 0;
            }
        }
    }
}

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <class T>
void print(const T& c) {
    for (auto x : c) std::cout << x << ' ';
}

std::vector<double> householder_method(const std::vector<double>& matrix_, int size) {
    std::vector<double> matrix = matrix_;
    for (int i = 0; i < size - 2; ++i) {
        clock_t start = clock();
        std::vector<double> vector_s(size - i - 2);
        std::vector<double> vector_u(size - i - 2);
        std::vector<double> vector_v(size);
        double gamma;
        for (int j = i + 2; j < size; ++j) {
            vector_s.at(j - i - 2) = matrix.at(j * size + i);
        }
        if (vector_norm2(vector_s) < EPS) {
            vector_v.at(i) = 1;
            gamma = 0.5;
        } else {
            double s_norm = vector_norm(vector_s);
            for (int j = 0; j < size - i - 2; ++j) {
                vector_u.at(j) = vector_s.at(j) / s_norm;
            }
        
            for (int j = 0; j < size; ++j) {
                if (j < i) {
                    vector_v.at(j) = 0;
                } else if(j > i) {
                    vector_v.at(j) = vector_u.at(j - i );
                } else if (abs(vector_u[0]) < EPS) {
                    vector_v.at(j) = 1;
                } else {
                    vector_v.at(j) =  sgn(vector_u[0]) * (1 + abs(vector_u[0]));
                }
                gamma = 1 + abs(vector_u[0]);
            }
        }

        // std::cout << "v: ";
        // print(vector_v);
        // std::cout << std::endl;


        // //init vector_u
        // for (int j = 0; j < size; ++j) {
        //     if (j < i + 2) {
        //         vector_u[j] = 0;
        //     } else {
        //         vector_u[j] = matrix[(j) * size + i];
        //     }
        // }
        // vector_u[i + 2] += vector_norm(vector_u);



        // print(matrix); std::cout << std::endl;
        std::vector<double> tmp_matrix(size * size);
        tmp_matrix = matrix_multiplication(vector_v, matrix, 1, size, size);
        //std::cout << "tmp_matrix_size_after_1_mul = " << tmp_matrix.size() << std::endl;
        tmp_matrix = matrix_multiplication(vector_v, tmp_matrix, size, 1, size);
        //std::cout << "tmp_matrix_size_after_2_mul = " << tmp_matrix.size() << std::endl;
        tmp_matrix = mul_matrix_by_number(tmp_matrix, 1 / gamma, size);
        //std::cout << "tmp_matrix_size_after_num_mul = " << tmp_matrix.size() << std::endl;
        matrix = matrix_subtraction(matrix, tmp_matrix, size);

        //matrix = matrix_subtraction(matrix, mul_matrix_by_number(matrix_multiplication(matrix_multiplication(matrix, vector_u, size, size, 1), vector_u, size, 1, size), 2/vector_norm(vector_u), size), size);
        // std::vector<double> householder_matrix = generate_household_matrix_from_vector_u(vector_u);

        // std::vector<double> right_down_block = get_right_down_block(matrix, size, size - 1 - i);
        // right_down_block = matrix_multiplication(householder_matrix, right_down_block, size - 1 - i, size - 1 - i, size - 1 - i);
        // right_down_block = matrix_multiplication(right_down_block, householder_matrix, size - 1 - i, size - 1 - i, size - 1 - i);
        // std::vector<double> left_down_block = get_left_down_block(matrix, size, i);
        // std::vector<double> left_top_block = get_left_top_block(matrix, size, i);
        // std::vector<double> right_top_block = get_left_down_block(matrix, size, i);
        // //left_down_block = matrix_multiplication(householder_matrix, left_down_block, size - 1 - i, size - 1 - i, 1 + i);
        // right_down_block = matrix_multiplication(right_top_block, householder_matrix, 1 + i, size - 1 - i, size - 1 - i);
        // //matrix = get_matrix_from_blocks(left_top_block, right_top_block, left_down_block, right_down_block, size, i);
    }
    //std::cout << matrix.size() << std::endl;
    find_the_zeros_of_the_matrix(matrix, size);
    //std::cout << "!\n";
    return matrix;
}