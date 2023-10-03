#include <iostream>
#include <cstdlib>
#include <algorithm> 
#include <ctime> 

//k,j,i

//i,k,j - fastest

int main() {
    clock_t start = clock();
    int n = 1024;
    double* a = new double(n*n);
    double* b = new double(n*n);
    double* c = new double(n*n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i * n + j] = 1;
            b[i * n + j] = 1;
            c[i * n + j] = 0;
        }
    }
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << seconds << std::endl;
}