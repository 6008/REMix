#include <stdlib.h>
#include <stdio.h>
#include "lapacke.h"
#include "const.h"
#include "tools.h"

int lse(double* a, double* b, lapack_int m, lapack_int n)
{
	int i = 0;
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
	/*
	LAPACKE_dgelsd_work(matrix_layout, m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork );
	*/
	LAPACKE_dgelsd(matrix_layout, m, n, nrhs, a, lda, b, ldb, s, rcond, rank);

	free(rank);
	free(s);
	free(work);
	free(iwork);
	return 0;
}

int initial_mu(double* intensity, int* samples, size_t sample_size, double* knots, size_t nknots, double* mu, size_t ncycles, size_t nreads)
{
	size_t read_length = ncycles * CHANNEL;
	size_t m = ncycles * sample_size;
	size_t n = nknots + 2;
	double* first = malloc(sizeof(double) * m);
	double* second = malloc(sizeof(double) * m);
	double* third = malloc(sizeof(double) * m);
	double* fourth = malloc(sizeof(double) * m);
	double* temp = malloc(sizeof(double) * CHANNEL);	
	double* a = malloc(sizeof(double) * m * n);
	double* a1 = malloc(sizeof(double) * m * n);
	double* a2 = malloc(sizeof(double) * m * n);
	double* a3 = malloc(sizeof(double) * m * n);
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	while (i < sample_size)
	{
		while (j < ncycles)
		{
			*(temp + 0) = *(intensity + (*(samples + i) * read_length + j * CHANNEL + 0));
			*(temp + 1) = *(intensity + (*(samples + i) * read_length + j * CHANNEL + 1));
			*(temp + 2) = *(intensity + (*(samples + i) * read_length + j * CHANNEL + 2));
			*(temp + 3) = *(intensity + (*(samples + i) * read_length + j * CHANNEL + 3));
			qsort(temp, CHANNEL, sizeof(double), cmp_double);
			*(first + i * ncycles + j) = *(temp + 0);
			*(second + i * ncycles + j) = *(temp + 1);
			*(third + i * ncycles + j) = *(temp + 2);
			*(fourth + i * ncycles + j) = *(temp + 3);
			j++;
		}
		j = 0;
		i++;
	}

	i = 0;
	j = 0;
	
	while (i < sample_size)
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
		while (i < sample_size)
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
	
	for (i = 0; i < m * n; i++)
	{
		*(a1 + i) = *(a + i);
		*(a2 + i) = *(a + i);
		*(a3 + i) = *(a + i);
	}
	
	i = 0;
	j = 0;
	k = 0;
	
	lse(a1, second, m, n);
	lse(a2, third, m, n);
	lse(a3, fourth, m, n);
	lse(a, first, m, n);
	
	while (k < n)
	{
		*(mu + k * CHANNEL * CHANNEL + 0) = *(first + k);
		*(mu + k * CHANNEL * CHANNEL + 1) = *(second + k);
		*(mu + k * CHANNEL * CHANNEL + 2) = *(third + k);
		*(mu + k * CHANNEL * CHANNEL + 3) = *(fourth + k);
		*(mu + k * CHANNEL * CHANNEL + 4) = *(second + k);
		*(mu + k * CHANNEL * CHANNEL + 5) = *(first + k);
		*(mu + k * CHANNEL * CHANNEL + 6) = *(fourth + k);
		*(mu + k * CHANNEL * CHANNEL + 7) = *(third + k);
		*(mu + k * CHANNEL * CHANNEL + 8) = *(third + k);
		*(mu + k * CHANNEL * CHANNEL + 9) = *(fourth + k);
		*(mu + k * CHANNEL * CHANNEL + 10) = *(first + k);
		*(mu + k * CHANNEL * CHANNEL + 11) = *(second + k);
		*(mu + k * CHANNEL * CHANNEL + 12) = *(fourth + k);
		*(mu + k * CHANNEL * CHANNEL + 13) = *(third + k);
		*(mu + k * CHANNEL * CHANNEL + 14) = *(second + k);
		*(mu + k * CHANNEL * CHANNEL + 15) = *(first + k);
		k++;
	}
	
	free(first);
	free(second);
	free(third);
	free(fourth);
	free(a);
	free(a1);
	free(a2);
	free(a3);
	free(temp);
	return 0;
}

int initial_sigma(size_t nknots, double** sigma)
{
	size_t i = 0;
	/*
	Because there are k_0, k_1, and other ks associated with nknots
	The total number of ks is (2 + nknots)
	*/
	size_t length = CHANNEL * (2 + nknots);
	size_t scale = length * length;

	while (i < CHANNEL)
	{
		size_t temp = 0;
		*(sigma + i) = malloc(sizeof(double) * scale);
		while (temp < scale)
		{
			*(*(sigma + i) + temp) = 0.0;
			temp++;
		}
		
		temp = 0;
		while (temp < length)
		{
			*(*(sigma + i) + temp * length + temp) = 1.0;
			temp++;
		}
		i++;
	}
	return 0;
}

int initial_pi(double* pi)
{
	size_t i = 0;
	while (i < CHANNEL)
	{
		*(pi + i) = 1.0 / CHANNEL;
		i++;
	}
	return 0;
}

int initial_v(double** v)
{
	size_t i = 0;
	size_t length = CHANNEL;
	size_t scale = length * length;

	while (i < CHANNEL)
	{
		size_t temp = 0;
		*(v + i) = malloc(sizeof(double) * scale);
		while (temp < scale)
		{
			*(*(v + i) + temp) = 0.0;
			temp++;
		}
		
		temp = 0;
		while (temp < length)
		{
			*(*(v + i) + temp * length + temp) = 1.0;
			temp++;
		}
		i++;
	}
	return 0;
}

int initial_probs(double* intensity, double* probs, size_t ncycles, size_t nreads)
{
	size_t i = 0;
	size_t j = 0;
	size_t total = ncycles * nreads * CHANNEL;
	size_t temp = 0;
	while (i < total)
	{
		*(probs + i) = 0.0;
		i++;
	}
	
	i = 0;
	
	while (i < nreads)
	{
		while (j < ncycles)
		{
			temp = i * ncycles * CHANNEL + j * CHANNEL;
			*(probs + temp + max_index(*(intensity + temp + 0), *(intensity + temp + 1), *(intensity + temp + 2), *(intensity + temp + 3))) = 1.0;
			j++;
		}
		j = 0;
		i++;
	}
	
	return 0;
}