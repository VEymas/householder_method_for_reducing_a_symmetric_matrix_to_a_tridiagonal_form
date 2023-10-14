#include <iostream>
#include <cstdlib>
#include <algorithm> 
#include <ctime> 


//i,k,j - fastest

void matmul(double * A, double * B, double * C, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i * n + j] = 0;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void block_matmul(double * A, double * B, double * C, int n, int block_size) {
    double * a = new double[block_size * block_size];
    double * b = new double[block_size * block_size];
    double * c = new double[block_size * block_size];
    for (int i = 0; i < (int) n / block_size; ++i) {
        for (int k = 0; k < (int) n / block_size; ++k) {
            for (int j = 0; j < (int) n / block_size; ++j) {
                for (int I = 0; I < block_size; ++I) {
                    for (int J = 0; J < block_size; ++J) {
                        a[I * block_size + J] = A[(block_size * i + I) * n + block_size * j + J];
                        b[I * block_size + J] = B[(block_size * j + I) * n + block_size * k + J];
                        c[I * block_size + J] = 0;
                    }
                }
                matmul(a, b, c, block_size);
                for (int I = 0; I < block_size; ++I) {
                    for (int J = 0; J < block_size; ++J) {
                        C[(block_size * i + I) * n + block_size * j + J] += c[I * block_size + J];
                    }
                }
                
            }
        }
    }

}

void print_matrix(double* matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << matrix[i * n + j] << "\t";
        }
        std::cout << std::endl;
    }
}

template <class T>
void print(const T& c) {
    for (auto x : c) std::cout << x << ' ';
}
#include <vector>
int main() {
    // int n = 10;
    // int block_size = 2;
    // double* a = new double[n * n];
    // double* b = new double[n * n];
    // double* c = new double[n * n];
    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         a[i * n + j] = (double)std::rand() / 10000000;;
    //         b[i * n + j] = (double)std::rand() / 10000000;;
    //     }
    // }
    // clock_t start = clock();
    // matmul(a, b, c, n);
    // clock_t end = clock();
    // double seconds_matmul = (double)(end - start) / CLOCKS_PER_SEC;
    // print_matrix(c, n);
    // std::cout << std::endl;
    // start = clock();
    // block_matmul(a, b, c, n, block_size);
    // end = clock();
    // double seconds_block_matmul = (double)(end - start) / CLOCKS_PER_SEC;
    // print_matrix(c, n);
    // std::cout << "matmul - " << seconds_matmul << "s block_matmul - " << seconds_block_matmul << "s" << std::endl;
    std::vector<double> test(10);
    print(test);
    return 0;
}