#ifndef BC_TOOLS_H
#define BC_TOOLS_H

/*
file format:
each line represents a read (4 * ncycles length)
each column represents a cycle

parameters:
file_name: the file contains intensity
intensity: record the intensity from file, should be well malloced
delima: delim between values
nreads: number of reads
ncycles: number of cycles

return:
0: success
1: something wrong
*/
int read_intensity(const char* file_name, double* intensity, const char delima, const size_t nreads, const size_t ncycles);

/*
simple quntile knots
begin should be less than end
nknots should be greater than 0
knots should be well malloced
*/
int simple_knots(double* knots, const int begin, const int end, const int nknots);


/*
sample N distinct numbers from begin(include) to end(exclude)
end - begin > N
*/
int sample(int* samples, int begin, int end, size_t nsample);

int read_samples(int* samples, size_t nsample, const char* file_name);


int cmp_double(const void* a, const void* b);

int max_index(double d1, double d2, double d3, double d4);

/*
Sum of d1 to d4 should be 1
*/
int rmultinom(double d1, double d2, double d3, double d4);

/*
begin include
end exclude
sum of begin to end -1 should be 1
*/
int rmultinom_array(double* arr, int begin, int end);

/*
which(Delta[,k,j]==value)
*/
int select_delta(int* delta, int* result, size_t nreads, size_t ncycles, size_t col_index, int value);
/*
matrix is m*n
the first element of index(_rows,_cols,) is the length of select array.
so the length of index(_rows,_cols,) is 1 + this value.
*/
int fetch_elems(double* matrix, size_t m, size_t n, int* index_rows, int* index_cols, double* result);
int fetch_cols(double* matrix, size_t m, size_t n, int* index, double* result);
int fetch_rows(double* matrix, size_t m, size_t n, int* index, double* result);

void print_matrix_rowmajor(char* desc, int m, int n, double* mat, int ldm);
void print_matrix_colmajor(char* desc, int m, int n, double* mat, int ldm);

void print_int_matrix_rowmajor(char* desc, int m, int n, int* mat, int ldm);

#endif