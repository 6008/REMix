#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "const.h"

char** str_split(char* str, const char delima)
{
	char** result = 0;
	size_t count = 0;
	char* tmp = str;
	char* last_comma = 0;
	char delim[2];
	delim[0] = delima;
	
	while (*tmp)
	{
		if(delima == *tmp)
		{
			count++;
			last_comma = tmp;
		}
		tmp++;
	}
	
	count += last_comma < (str + strlen(str) - 1);
	count ++;
	
	result = malloc(sizeof(char *) * count);

	if (result)
    {
        size_t idx  = 0;
        char* token = strtok(str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

int read_intensity(const char* file_name, double* intensity, const char delima, const size_t nreads, const size_t ncycles)
{
	FILE* fp = fopen(file_name, "r");
	char* line = 0;
	size_t len = 0;
	size_t nlines = 0;
	ssize_t read;
	
	if (fp == NULL)
		return 1;
	
	while ((read = getline(&line, &len, fp)) != -1)
	{
		char** values = str_split(line, delima);
		size_t idx = 0;
		while (idx < ncycles * CHANNEL)
		{
			*(intensity + nlines * ncycles * CHANNEL + idx) = atof(*(values + idx));
			idx++;
		}
		nlines++;
	}
	
	fclose(fp);
	
	return 0;
}

int simple_knots(double* knots, const int begin, const int end, const int nknots)
{
	size_t i = 0;
	if ((nknots < 1) || (begin >= end))
		return 1;
	while (i < nknots)
	{
		*(knots + i) = (end - begin + 0.0) * (i + 1) / (nknots + 1) + begin;
		i++;
	}
	return 0;
}

int sample(int* samples, int begin, int end, size_t nsample)
{
	size_t i = 0;
	if ((nsample <= 0) || (begin >= end) || (nsample > end - begin))
		return 1;
	while (i < nsample)
	{
		*(samples + i) = (rand() % (end - begin)) + begin;
		i++;
	}
	return 0;
}

int read_samples(int* samples, size_t nsample, const char* file_name)
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
		*(samples + i) = atoi(line) - 1;
		i++;
	}
	fclose(fp);
	return 0;
}

int cmp_double(const void* a, const void* b)
{
	if (*(double*)b > *(double*)a)
		return 1;
	else
		return -1;
}

int max_index(double d1, double d2, double d3, double d4)
{
	if ((d1 > d2) && (d1 > d3) && (d1 > d4))
		return 0;
	if ((d2 > d3) && (d2 > d4) && (d2 > d1))
		return 1;
	if ((d3 > d4) && (d3 > d1) && (d3 > d2))
		return 2;
	return 3;
}

int rmultinom(double d1, double d2, double d3, double d4)
{
	double r = (double)rand() / (double)RAND_MAX;
	if (r > d1 + d2 + d3)
		return 3;
	if (r > d1 + d2)
		return 2;
	if (r > d1)
		return 1;
	return 0;
}

int rmultinom_array(double* arr, int begin, int end)
{
	double r = (double)rand() / (double)RAND_MAX;
	double* temp = malloc(sizeof(double) * (end - begin));
	int i = begin;
	int result = begin;
	*(temp + 0) = 0.0;
	for (i = begin + 1; i < end; i++)
	{
		*(temp + i -  begin) = *(temp + i - begin - 1) + *(arr + i);
	}
	for ( i = end - 1; i > begin; i++)
	{
		if (r > *(temp + i))
		{
			result = i;
		}
	}
	free(temp);
	return result;
} 

int select_delta(int* delta, int* result, size_t nreads, size_t ncycles, size_t col_index, int value)
{
	size_t i = 0;
	size_t j = 1;
	
	for ( i = 0; i < nreads; i++)
	{
		if (*(delta + i * ncycles * CHANNEL + col_index) == value)
		{
			*(result + j) = i;
			j++;
		}
	}
	*(result + 0) = j - 1;
	return 0;
}

int fetch_rows(double* matrix, size_t m, size_t n, int* index, double* result)
{
	size_t i;
	size_t j;
	for (i = 0; i < *(index + 0); i++)
	{
		for (j = 0; j < n; j++)
		{
			*(result + i * n + j) = *(matrix + *(index + i + 1) * n + j);
		}
	}
	return 0;
}

int fetch_cols(double* matrix, size_t m, size_t n, int* index, double* result)
{
	size_t i;
	size_t j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < *(index + 0); j++)
		{
			*(result + i * (*(index + 0)) + j) = *(matrix + i * n + *(index + j + 1));
		}
	}
	return 0;
}

int fetch_elems(double* matrix, size_t m, size_t n, int* index_rows, int* index_cols, double* result)
{
	size_t i;
	size_t j;
	for (i = 0; i < *(index_rows + 0); i++)
	{
		for (j = 0; j < *(index_cols + 0); j++)
		{
			*(result + i * (*(index_cols + 0)) + j) = *(matrix + *(index_rows + i + 1) * n + *(index_cols + j + 1));
		}
	}
	return 0;
}


void print_matrix_rowmajor(char* desc, int m, int n, double* mat, int ldm)
{
	int i, j;
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

void print_matrix_colmajor(char* desc, int m, int n, double* mat, int ldm)
{
	int i, j;
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

void print_int_matrix_rowmajor(char* desc, int m, int n, int* mat, int ldm)
{
	int i, j;
	printf("\n %s\n", desc);
	
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf( " %8d", mat[i*ldm+j] );
		}
		printf( "\n" );
	}
}