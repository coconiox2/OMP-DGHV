#ifndef __MATRIX_H_
#define __MATRIX_H_

int get_matrix_int(int  **output, int size_n, int min_val, int max_val);
int print_matrix_int(int *matrix, int size_n);
int matrix_mul_and_compare(int *matrix1, int *matrix2, int size, int *compare, float *duration);

#endif // __MATRIX_H_