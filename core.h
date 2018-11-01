








/*
This is the only statement in this file.
*/
int mcem(double* intensity, double* zdenom, double* mu, double** sigma, double** v, double* pi, double* probs, double* knots, size_t nreads, size_t ncycles, size_t niter, size_t nknots, double q);





int getdenom(double* intensity, double* denom, double* bigx, double* mu, double* pi, double** v, double** sigma, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int update_v(double** egamma, double** covgamma, double** v, double* intensity, double* bigx, int* delta, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int update_sigma(double** sigma, double** covgamma, double** egamma, double* mu, size_t nreads, size_t bigx_cols);
int update_mu(double** egamma, double* mu, size_t nreads, size_t bigx_cols);
int update_pi(int* delta, double* pi, size_t nreads, size_t ncycles);
int posterior_means(int** index_group, double** egamma, double** covgamma, double** v, double** sigma, double* intensity, double* bigx, double* mu, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int init_index_group(int** index_group, int* delta, size_t nreads, size_t ncycles);
int bigx_row_major(double* bigx, size_t ncycles, size_t nknots, double* knots);
int update_delta(double* denom, int* delta, double q, size_t nreads, size_t ncycles);
