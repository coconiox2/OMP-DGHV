#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define CL_DEBUG	0

int get_matrix_int(int  **output, int size_n, int min_val, int max_val)
{
	int ret = 0;
	*output = (int *)calloc(size_n * size_n, sizeof(int));
	srand(time(NULL));
	int i, j;
	int range = max_val - min_val + 1;

	for (i = 0; i < size_n; i++)
	{
		for (j = 0; j < size_n; j++)
		{
			(*output)[i * size_n + j] = (rand() % range) + min_val; 
		}
	}

	return ret;
}

int print_matrix_int(int *matrix, int size_n)
{
	int ret = 0;
	int i, j;
	for (i = 0; i < size_n; i++)
	{
		for (j = 0; j < size_n; j++)
		{
			printf("%-4d ", matrix[i * size_n + j]);
		}
		printf("\n");
	}
	printf("\n");
	getchar();

	return ret;
}

int matrix_mul_and_compare(int *matrix1, int *matrix2, int size, int *compare, float *duration)
{
	int i, j, k, element = 0;
	*duration = 0;
	clock_t clocks = -clock();

	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				element += matrix1[i * size + k] * matrix1[k * size + j];
			}
#if CL_DEBUG
			if (element != compare[i * size + j]) 
			{
				fprintf(stdout, "FAILURE!\n");
				return -1;
			}
			else
			{
				element = 0;
			}
#else
			compare[i * size + j] = element;
			element = 0;
#endif
		}
	}

	clocks += clock();
	*duration = ((float)clocks) / CLOCKS_PER_SEC;

	fprintf(stderr, "SUCCESS!\n");

	return 0;
}