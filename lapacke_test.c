#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

int main_test_dgetrf();
int main_test_dgemm_special();
int main_test_simple_lse();
int main_test_lse(char* file_name);
int main_test_lapack();
int main_test_dgemm();
int main_test_free();
int main_test_inverse();
void print_matrix_rowmajor(char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm);
void print_matrix_colmajor(char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm);

int lse(double* a, double* b, lapack_int m, lapack_int n)
{
	int matrix_layout = LAPACK_COL_MAJOR;
	int lwork = m * n;
	double rcond = 0.0001;
	lapack_int nrhs = 1;
	lapack_int lda = m;
	lapack_int ldb = m;
	lapack_int* rank = malloc(sizeof(lapack_int) * 1);
	lapack_int* iwork = malloc(sizeof(lapack_int) * lwork);
	double* work = malloc(sizeof(double) * lwork);
	double* s = malloc(sizeof(double) * n);
	LAPACKE_dgelsd_work(matrix_layout, m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork );
	free(rank);
	free(s);
	free(work);
	free(iwork);
	return 0;
}

int read_samples(double* samples, size_t nsample, const char* file_name)
{
	FILE* fp = fopen(file_name, "r");
	char* line = 0;
	size_t len = 0;
	size_t i = 0;
	ssize_t read;
	
	if (fp == NULL)
		return 1;
	
	while ((read = getline(&line, &len, fp)) != -1)
	{
		*(samples + i) = atof(line);
		i++;
	}
	fclose(fp);
	return 0;
}

int main(int argc, char** argv)
{
	printf("%f\n", pow(2,-0.5));
	return main_test_dgetrf();
}

/*
lapack_int LAPACKE_dgetri( int matrix_layout, lapack_int n, double* a,
                           lapack_int lda, const lapack_int* ipiv );
*/

int main_test_inverse()
{
	int matrix_layout = LAPACK_ROW_MAJOR;
	double matrix[6][6] = {{42,3,4,5,6,7}, {5,24,4,3,2,1}, {3,7,32,1,2,3},
						   {0,0,0,18,0,0}, {0,0,0,0,31,0}, {7,6,5,4,3,52}};
	int* ipiv = malloc(sizeof(int) * 6);
	print_matrix_rowmajor("Result", 6, 6, *matrix, 6);
	LAPACKE_dgetrf(matrix_layout, 6, 6, *matrix, 6, ipiv);
	print_matrix_rowmajor("Result", 6, 6, *matrix, 6);
	LAPACKE_dgetri(matrix_layout, 6, *matrix, 6, ipiv);
	print_matrix_rowmajor("Result", 6, 6, *matrix, 6);
	return 0;
}

int main_test_dgemm_special()
{
	int m = 4;
	int k = 3;
	int n = 2;
	double alpha = 1;
	double beta = 1;
	double a[20] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5};
	double b[12] = {1,1,1,2,2,2,3,3,3,4,4,4};
	double c[15] = {1,1,1,2,2,2,3,3,3,4,4,4,5,5,5};
	print_matrix_rowmajor("A", 5, 4, a, 4);
	print_matrix_rowmajor("B", 4, 3, b, 3);
	print_matrix_rowmajor("C", 5, 3, c, 3);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	m, n, k, alpha, a + 5, 4, b, 3, beta, c, 3);
	print_matrix_rowmajor("C", 5, 3, c, 3);
	return 0;
}

/*
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
*/

int main_test_dgemm()
{
	int m = 4;
	int k = 3;
	int n = 2;
	double alpha = 1;
	double beta = 0;
	double* a = malloc(sizeof(double) * m * k);
	double* b = malloc(sizeof(double) * n * k);
	double* c = malloc(sizeof(double) * m * n);
	int i = 0;
	int j = 0;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			*(a + i * k + j) = 1;
			printf("%f\t", *(a + i * k + j));
		}
		printf("\n");
	}
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(b + i * n + j) = i * n + j;
			printf("%f\t", *(b + i * n + j));
		}
		printf("\n");
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(c + i * n + j) = 0;
		}
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, k, b, n, beta, c, n);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%f\t", *(c + i * n + j));
		}
		printf("\n");
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(c + i * n + j) = 0;
		}
	}
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, a, k, b, k, beta, c, n);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%f\t", *(c + i * n + j));
		}
		printf("\n");
	}
	return 0;
}

int main_test_free()
{
	double* test = malloc(sizeof(double) * 5);
	double* test2 = malloc(sizeof(double) * 5);
	size_t i = 0;
	for (i = 0; i < 5; i++)
	{
		*(test + i) = 0;
		*(test2 + i) = i;
	}
	free(test);
	test = test2;
	for (i = 0; i < 5; i++)
	{
		printf("%f\t", *(test + i));
	}
	free(test);
	return 0;
}

int main_test_simple_lse()
{
	double A[20] = {1,1,2,3,1,4,5,1,1,3,5,2,1,4,1,4,1,2,5,3};
	double b[5] = {-10, 14, 18, 14, 16};
	lapack_int info,m,n,lda,ldb,nrhs;
	m = 5;
	n = 4;
	nrhs = 1;
	lda = 4;
	ldb = 1;
	
	print_matrix_rowmajor("A", m, n, A, lda);
	print_matrix_rowmajor("b", m, nrhs, b, ldb);
	info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,A,lda,b,ldb);
	print_matrix_rowmajor("b", m, nrhs, b, ldb);
	return 0;
}

int main_test_lse(char* file_name)
{
	/*
lapack_int LAPACKE_dgelsd_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* s,
                                double rcond, lapack_int* rank, double* work,
                                lapack_int lwork, lapack_int* iwork );
	*/
	int matrix_layout = LAPACK_COL_MAJOR;
	int nknots = 5;
	int ncycles = 101;
	int nreads = 100;
	int lwork = ncycles * nreads * 10;
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	lapack_int m = ncycles * nreads;
	lapack_int n = 2 + nknots;
	lapack_int nrhs = 1;
	lapack_int lda = m;
	lapack_int ldb = m;
	lapack_int* rank = malloc(sizeof(lapack_int) * 1);
	lapack_int* iwork = malloc(sizeof(lapack_int) * lwork);
	double rcond = 0.0001;
	double* knots = malloc(sizeof(double) * nknots);
	double* b = malloc(sizeof(double) * m);
	double* a = malloc(sizeof(double) * n * m);
	double* s = malloc(sizeof(double) * n);
	double* work = malloc(sizeof(double) * lwork);
	*(knots + 0) = 17.5;
	*(knots + 1) = 34;
	*(knots + 2) = 50.5;
	*(knots + 3) = 67;
	*(knots + 4) = 83.5;
	read_samples(b, m, file_name);
	i = 0;
	while (i < nreads)
	{
		while(j < ncycles)
		{
			*(a + (i * ncycles + j)) = 1;
			*(a + (m + i * ncycles + j)) = j;
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
				*(a + (m * 2 + k * m + i * ncycles + j)) = (j > *(knots + k)) ? j - *(knots + k) : 0;
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
	/*
	print_matrix_colmajor("A", m, n, a, lda);

	LAPACKE_dgelsd_work(matrix_layout, m, n, nrhs, a, lda, b, ldb, s, rcond, rank,  work, lwork, iwork );

	LAPACKE_dgels_work(matrix_layout, 'N', m, n, nrhs, a, lda, b, ldb, work, lwork);
	*/
	lse(a, b, m, n);
	while ( k < n )
	{
		printf("%f\n", *(b + k));
		k++;
	}
	
	free(a);
	free(b);
	free(rank);
	free(s);
	free(work);
	free(iwork);
	free(knots);
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
	int nreads = 50;
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
	double* x = malloc(sizeof(double) * (m + 1));
	double* d = malloc(sizeof(double) * n);
	double* b = malloc(sizeof(double) * n * p);
	double* a = malloc(sizeof(double) * n * m);
	
	*(knots + 0) = 1.75;
	*(knots + 1) = 3.4;
	*(knots + 2) = 5.05;
	*(knots + 3) = 6.7;
	*(knots + 4) = 8.35;
	
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
		*(b + (i * p + i)) = 1;
		i++;
	}
	i = 0;
	j = 0;
	
	while (i < nreads)
	{
		while(j < ncycles)
		{
			*(a + (i * ncycles + j)) = 1;
			*(a + (n + i * ncycles + j)) = j;
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
				*(a + (n * 2 + k * n + i * ncycles + j)) = (j > *(knots + k)) ? j - *(knots + k) : 0;
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
/*
	while (k < m)
	{
		while (i < n)
		{
			printf("%f\t", *(a + k * n + i));
			i++;
		}
		i = 0;
		printf("\n");
		k++;
	}
	k = 0;
	
*/
	
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
	return 0;
}

void print_matrix_rowmajor(char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm)
{
	lapack_int i, j;
	printf("\n %s\n", desc);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf( " %6.2f", mat[i*ldm+j] );
		}
		printf( "\n" );
	}
}

void print_matrix_colmajor(char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm)
{
	lapack_int i, j;
	printf("\n %s\n", desc);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf( " %6.2f", mat[i+ldm*j] );
		}
		printf( "\n" );
	}
}

/*
lapack_int LAPACKE_dgetrf2( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );
*/

int main_test_dgetrf()
{
	int m = 5;
	int n = 5;
	int lda = 5;
	double a[5][5] = {{1,3,5,7,9},{5,4,3,2,1},{7,3,4,6,8},{2,9,5,9,9},{3,5,7,6,1}};
	int *ipiv = malloc(sizeof(int) * 5);
	print_matrix_rowmajor("Before", m, n, *a, lda);
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, *a, lda, ipiv);
	print_matrix_rowmajor("After", m, n, *a, lda);
	printf("%d\t%d\t%d\t%d\t%d\n", *(ipiv + 0), *(ipiv + 1),*(ipiv + 2),*(ipiv + 3),*(ipiv + 4));
	printf("%f\n", a[0][0] * a[1][1] * a[2][2] * a[3][3] * a[4][4]);
}