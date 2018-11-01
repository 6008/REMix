#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "tools.h"
#include "core.h"

int main_test_input();
int main_test_initial();
int main_test_core();

int main()
{
	srand(time(NULL));
	return main_test_core();
}

int main_test_core()
{
	char* file_name = "test.txt";
	char* prob_name = "probs.txt";
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t ncycles = 101;
	size_t nreads = 100;
	size_t scale = CHANNEL * ncycles * nreads;
	size_t nknots = 5;
	size_t niter = 10;
	double q = 0.01;
	int begin = 1;
	int end = ncycles - 1;
	int sample_size = nreads;
	int* samples = malloc(sizeof(int) * sample_size);
	double** sigma = malloc(sizeof(double*) * CHANNEL);
	double** v = malloc(sizeof(double*) * CHANNEL);
	double* pi = malloc(sizeof(double) * CHANNEL);
	double* knots = malloc(sizeof(double) * nknots);
	double* intensity = malloc(sizeof(double) * scale);
	double* mu = malloc(sizeof(double) * CHANNEL * CHANNEL * (2 + nknots));
	double* zdenom = malloc(sizeof(double) * scale);
	double* probs = malloc(sizeof(double) * scale);
	read_intensity(file_name, intensity, ' ', nreads, ncycles);
	read_intensity(prob_name, probs, ' ', nreads, ncycles);
	simple_knots(knots, begin, end, nknots);
	for (i = 0; i < sample_size; i++)
	{
		*(samples + i) = i;
	}
	initial_mu(intensity, samples, sample_size, knots, nknots, mu, ncycles, nreads);
	initial_sigma(nknots, sigma);
	initial_pi(pi);
	initial_v(v);
	/*
	initial_probs(intensity, probs, ncycles, nreads);
	*/
	mcem(intensity, zdenom, mu, sigma, v, pi, probs, knots, nreads, ncycles, niter, nknots, q);
	
	for (k = 0; k < CHANNEL; k++)
	{
		free(*(sigma + k));
		free(*(v + k));
	}
	free(samples);
	free(sigma);
	free(v);
	free(pi);
	free(knots);
	free(intensity);
	free(mu);
	free(zdenom);
	free(probs);
	return 0;
}

int main_test_initial()
{
	char* file_name = "test.txt";
	char* sample_file = "samples.txt";
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t ncycles = 101;
	size_t nreads = 100;
	size_t scale = CHANNEL * ncycles * nreads;
	size_t nknots = 5;
	int begin = 1;
	int end = ncycles - 1;
	int sample_size = nreads;
	int* samples = malloc(sizeof(int) * sample_size);
	double* knots = malloc(sizeof(double) * nknots);
	double* intensity = malloc(sizeof(double) * scale);
	double* mu = malloc(sizeof(double) * CHANNEL * CHANNEL * (2 + nknots));
	read_intensity(file_name, intensity, ' ', nreads, ncycles);
	simple_knots(knots, begin, end, nknots);
	/*
	read_samples(samples, sample_size, sample_file);
	*/
	for (i = 0; i < sample_size; i++)
	{
		*(samples + i) = i;
	}
	initial_mu(intensity, samples, sample_size, knots, nknots, mu, ncycles, nreads);
	printf("Print knots:\n");
	while (i < nknots)
	{
		printf("%f\t", *(knots + i));
		i++;
	}
	printf("\n");
	i = 0;
	printf("Print Mu:\n");
	while (i < CHANNEL * (2 + nknots))
	{
		while (j < CHANNEL)
		{
			printf("%f\t", *(mu + i * CHANNEL + j));
			j++;
		}
		j = 0;
		printf("\n");
		i++;
	}
	free(mu);
	free(intensity);
	free(samples);
	free(knots);
	return 0;
}

int main_test_input()
{
	char* file_name = "data.txt";
	size_t ncycles = 1;
	size_t nreads = 3;
	size_t scale = 4 * ncycles * nreads;
	double* intensity = malloc(sizeof(double) * scale);
	char delima = '\t';
	read_intensity(file_name, intensity, delima, nreads, ncycles);
	
	while (scale > 0)
	{
		printf("%f\n", *(intensity + scale - 1));
		scale--;
	}
	
	free(intensity);
	return 0;
}

