#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "tools.h"
#include "initial.h"
#include "lapacke.h"

int main_test_lapack();
int main_test_qsort();
int main_test_split();
int main_test_getline();
int main_test_knots();
int main_test_initial();
int main_test_fetch();

int main()
{
	return main_test_fetch();
}

int main_test_fetch()
{
	double matrix[4][5] = {{0,1,2,3,4}, 
						   {5,6,7,8,9},
						   {1,3,5,7,9},
						   {2,4,6,8,0}};
	int m = 4;
	int n = 5;
	int index_row[3] = {2, 1, 3};
	int index_col[4] = {3, 0, 2, 4};
	double* result = malloc(sizeof(double) * 6);
	double* result_row = malloc(sizeof(double) * 10);
	double* result_col = malloc(sizeof(double) * 12);
	fetch_elems(*matrix, m, n, index_row, index_col, result);
	fetch_cols(*matrix, m, n, index_col, result_col);
	fetch_rows(*matrix, m, n, index_row, result_row);
	print_matrix_rowmajor("result", 2, 3, result, 3);
	print_matrix_rowmajor("result_row", 2, 5, result_row, 5);
	print_matrix_rowmajor("result_col", 4, 3, result_col, 3);
	return 0;
}

int main_test_lapack()
{
	/*
	lapack_int LAPACKE_dggglm( int matrix_layout, lapack_int n, lapack_int m,
							lapack_int p, double* a, lapack_int lda, double* b,
							lapack_int ldb, double* d, double* x, double* y );
	*/
	int matrix_layout = LAPACK_COL_MAJOR;
	int nknots = 5;
	int ncycles = 101;
	int nreads = 100;
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	lapack_int n = ncycles * nreads;
	lapack_int m = 2 + nknots;
	lapack_int p = n;
	lapack_int lda = n;
	lapack_int ldb = n;
	double* knots = malloc(sizeof(double) * nknots);
	double* y = malloc(sizeof(double) * p);
	double* x = malloc(sizeof(double) * m);
	double* d = malloc(sizeof(double) * n);
	double* b = malloc(sizeof(double) * n * p);
	double* a = malloc(sizeof(double) * n * m);
	
	*(knots + 0) = 17.5;
	*(knots + 1) = 34;
	*(knots + 2) = 50.5;
	*(knots + 3) = 67;
	*(knots + 4) = 83.5;
	
	while (i < n)
	{
		*(d + i) = 1;
		i++;
	}
	i = 0;
	while (i < n * p)
	{
		*(b + i) = 0;
		i++;
	}
	i = 0;
	while (i < n)
	{
		while (j < p)
		{
			*(b + (i * p + j)) = 1;
			j++;
		}
		j = 0;
		i++;
	}
	i = 0;
	j = 0;
	
	while (i < nreads)
	{
		while(j < ncycles)
		{
			*(a + (i * ncycles + j)) = j;
			j++;
		}
		j = 0;
		i++;
	}
	i = 0;
	j = 0;
	while (k < nknots)
	{
		while (i < nreads)
		{
			while(j < ncycles)
			{
				*(a + (n + k * n + i * ncycles + j)) = (j > *(knots + k)) ? j - *(knots + k) : 0;
				j++;
			}
			j = 0;
			i++;
		}
		i = 0;
		k++;
	}
	i = 0;
	j = 0;
	k = 0;
	
	LAPACKE_dggglm(matrix_layout, n, m, p, a, lda, b, ldb, d, x, y );
	
	while ( k < m )
	{
		printf("%f\n", *(x + k));
		k++;
	}
	
	free(a);
	free(b);
	free(d);
	free(x);
	free(y);
	free(knots);
}

int main_test_qsort()
{
	double* test = malloc(sizeof(double) * CHANNEL);
	*(test + 0) = 3.4;
	*(test + 1) = 2.7;
	*(test + 2) = 5.6;
	*(test + 3) = 4.2;
	qsort(test, CHANNEL, sizeof(double), cmp_double);
	printf("%f\t%f\t%f\t%f\n", *(test + 0), *(test + 1), *(test + 2), *(test + 3));
}

int main_test_initial()
{
	double** sigma = malloc(sizeof(double*) * CHANNEL);
	double** v = malloc(sizeof(double*) * CHANNEL);
	double* pi = malloc(sizeof(double) * CHANNEL);
	size_t nknots = 1;
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	initial_sigma(nknots, sigma);
	initial_pi(pi);
	initial_v(v);
	printf("Print Sigma:\n");
	while (i < CHANNEL)
	{
		j = 0;
		while (j < CHANNEL * (nknots + 2))
		{
			k = 0;
			while (k < CHANNEL * (nknots + 2))
			{
				printf("%f\t", *(*(sigma + i) + CHANNEL * (nknots + 2) * j + k));
				k++;
			}
			printf("\n");
			j++;
		}
		printf("\n\n");
		free(*(sigma + i));
		i++;
	}
	free(sigma);
	printf("Print Pi:\n");
	i = 0;
	j = 0;
	k = 0;
	while (i < CHANNEL)
	{
		printf("%f\t", *(pi + i));
		i++;
	}
	printf("\n");
	free(pi);
	i = 0;
	printf("Print V:\n");
	while (i < CHANNEL)
	{
		j = 0;
		while (j < CHANNEL)
		{
			k = 0;
			while (k < CHANNEL)
			{
				printf("%f\t", *(*(v + i) + CHANNEL * j + k));
				k++;
			}
			printf("\n");
			j++;
		}
		printf("\n\n");
		free(*(v + i));
		i++;
	}
	free(v);
	return 0;
}

int main_test_knots()
{
	int begin = 1;
	int end = 100;
	int nknots = 5;
	double* knots = malloc(sizeof(double) * nknots);
	simple_knots(knots, begin, end, nknots);
	while (nknots > 0)
	{
		printf("%f\t", *(knots + nknots - 1));
		nknots--;
	}
	printf("\n");
	if(knots)
		free(knots);
	return 0;
}

int main_test_split()
{
    char months[] = "JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC";
    char** tokens;

    printf("months=[%s]\n\n", months);

    tokens = str_split(months, ',');

    if (tokens)
    {
        int i;
        for (i = 0; *(tokens + i); i++)
        {
            printf("month=[%s]\n", *(tokens + i));
            free(*(tokens + i));
        }
        printf("\n");
        free(tokens);
    }

    return 0;
}

int main_test_getline()
{
	FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen("data.txt", "r");
    if (fp == NULL)
        return 1;

    while ((read = getline(&line, &len, fp)) != -1) {
        printf("Retrieved line of length %zu :\n", read);
        printf("%d\t%s", len, line);
    }

    fclose(fp);
    if (line)
        free(line);
	return 0;
}