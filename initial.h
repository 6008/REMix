#ifndef BC_INITIAL_H
#define BC_INITIAL_H

/*
Get initial values of Mu
*/
int initial_mu(double* intensity, int* samples, size_t sample_size, double* knots, size_t nknots, double* mu, size_t ncycles, size_t nreads);

int initial_sigma(size_t nknots, double** sigma);

int initial_pi(double* pi);

int initial_v(double** v);

#endif