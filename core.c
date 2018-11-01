#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include "tools.h"
#include "const.h"

int getdenom(double* intensity, double* denom, double* bigx, double* mu, double* pi, double** v, double** sigma, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int update_v(double** egamma, double** covgamma, double** v, double* intensity, double* bigx, int* delta, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int update_sigma(double** sigma, double** covgamma, double** egamma, double* mu, size_t nreads, size_t bigx_cols);
int update_mu(double** egamma, double* mu, size_t nreads, size_t bigx_cols);
int update_pi(int* delta, double* pi, size_t nreads, size_t ncycles);
int posterior_means(int** index_group, double** egamma, double** covgamma, double** v, double** sigma, double* intensity, double* bigx, double* mu, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols);
int init_index_group(int** index_group, int* delta, size_t nreads, size_t ncycles);
int bigx_row_major(double* bigx, size_t ncycles, size_t nknots, double* knots);
int update_delta(double* denom, int* delta, double q, size_t nreads, size_t ncycles);
double det(double* matrix, size_t length, int* ipiv);

int mcem(double* intensity, double* zdenom, double* mu, double** sigma, double** v, double* pi, double* probs, double* knots, size_t nreads, size_t ncycles, size_t niter, size_t nknots, double q)
{
	size_t keep_it = niter - niter / 2;
	size_t it = 0;
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t bigx_rows = CHANNEL * ncycles;
	size_t bigx_cols = CHANNEL * (nknots + 2);
	size_t bigx_length = bigx_rows * bigx_cols;
	size_t sigma_length = bigx_cols * bigx_cols;
	size_t v_length = CHANNEL * CHANNEL;
	size_t int_length = bigx_rows * nreads;
	double* kfn = malloc(sizeof(double) * nknots * ncycles);
	double* bigx = malloc(sizeof(double) * bigx_rows * bigx_cols); 
	double* pi_avg = malloc(sizeof(double) * CHANNEL);
	double* mu_avg = malloc(sizeof(double) * bigx_cols * CHANNEL);
	double** sigma_avg = malloc(sizeof(double*) * CHANNEL);
	double** v_avg = malloc(sizeof(double*) * CHANNEL);
	double** egamma = malloc(sizeof(double*) * CHANNEL);
	double** covgamma = malloc(sizeof(double*) * CHANNEL * nreads);
	/*Same structure with intensity*/
	int* delta = malloc(sizeof(int) * int_length);
	int** index_group = malloc(sizeof(int*) * CHANNEL * nreads);
	
	for (i = 0; i < CHANNEL * nreads; i++)
	{
		*(covgamma + i) = malloc(sizeof(double) * bigx_cols * bigx_cols);
	}
	
	for (i = 0; i < CHANNEL; i++)
	{
		*(egamma + i) = malloc(sizeof(double) * nreads * bigx_cols);
		*(sigma_avg + i) = malloc(sizeof(double) * sigma_length);
		*(v_avg + i) = malloc(sizeof(double) * v_length);
	}
	i = 0;
	
	
	/*
	Init bigx
	*/
	for (i = 0; i < bigx_length; i++)
	{
		*(bigx + i) = 0.0;	
	}

	/*bigx_col_major(bigx, ncycles, nknots, knots);*/
	bigx_row_major(bigx, ncycles, nknots, knots);
	
	/*
	Init delta
	*/
	for (i = 0; i < int_length; i++)
	{
		*(delta + i) = (*(probs + i) > 0.5) ? 1 : 0;
	}
	
	/*
	Getting the indices
	*/
	init_index_group(index_group, delta, nreads, ncycles);
	
	/*
	
	for (i = 0; i < nreads; i++)
	{
		for (k = 0; k < CHANNEL; k++)
		{
			print_int_matrix_rowmajor("index_group", 1, *(*(index_group + i * CHANNEL + k) + 0) + 1, *(index_group + i * CHANNEL + k), *(*(index_group + i * CHANNEL + k) + 0) + 1);
		}
	}
	*/
	
	printf("pmeans\n");
	
	posterior_means(index_group, egamma, covgamma, v, sigma, intensity, bigx, mu, nreads, ncycles, bigx_rows, bigx_cols);

	
	printf("updating pi\n");
	update_pi(delta, pi, nreads, ncycles);
	printf("updating mu\n");
	update_mu(egamma, mu, nreads, bigx_cols);
	printf("updating sigma\n");
	update_sigma(sigma, covgamma, egamma, mu, nreads, bigx_cols);
	printf("updating v\n");
	update_v(egamma, covgamma, v, intensity, bigx, delta, nreads, ncycles, bigx_rows, bigx_cols);
	for (i = 0; i < CHANNEL; i++)
	{
		print_matrix_rowmajor("v", CHANNEL, CHANNEL, *(v + i), CHANNEL);
	}
	 
	printf("enter iterations\n"); 
	 
	for (it = 0; it < niter / 2; it++)
	{
		getdenom(intensity, zdenom, bigx, mu, pi, v, sigma, nreads, ncycles, bigx_rows, bigx_cols);
		update_delta(zdenom, delta, q, nreads, ncycles);
		init_index_group(index_group, delta, nreads, ncycles);
		posterior_means(index_group, egamma, covgamma, v, sigma, intensity, bigx, mu, nreads, ncycles, bigx_rows, bigx_cols);
		update_pi(delta, pi, nreads, ncycles);
		update_mu(egamma, mu, nreads, bigx_cols);
		update_sigma(sigma, covgamma, egamma, mu, nreads, bigx_cols);
		update_v(egamma, covgamma, v, intensity, bigx, delta, nreads, ncycles, bigx_rows, bigx_cols);		
	}
		
	for (i = 0; i < CHANNEL; i++)
	{
		*(pi_avg + i) = 0;
	}
	for (i = 0; i < bigx_cols * CHANNEL; i++)
	{
		*(mu_avg + i) = 0;
	}
	for (k = 0; k < CHANNEL; k++)
	{
		for (i = 0; i < sigma_length; i++)
		{
			*(*(sigma_avg + k) + i) = 0;
		}
		for (i = 0; i < v_length; i++)
		{
			*(*(v_avg + k) + i) = 0;
		}
	}	
	

	for (it = niter / 2; it < niter; it++)
	{
		getdenom(intensity, zdenom, bigx, mu, pi, v, sigma, nreads, ncycles, bigx_rows, bigx_cols);
		update_delta(zdenom, delta, q, nreads, ncycles);
		init_index_group(index_group, delta, nreads, ncycles);
		posterior_means(index_group, egamma, covgamma, v, sigma, intensity, bigx, mu, nreads, ncycles, bigx_rows, bigx_cols);
		update_pi(delta, pi, nreads, ncycles);
		update_mu(egamma, mu, nreads, bigx_cols);
		update_sigma(sigma, covgamma, egamma, mu, nreads, bigx_cols);
		update_v(egamma, covgamma, v, intensity, bigx, delta, nreads, ncycles, bigx_rows, bigx_cols);		
		for (i = 0; i < CHANNEL; i++)
		{
			*(pi_avg + i) += *(pi + i);
		}
		for (i = 0; i < bigx_cols * CHANNEL; i++)
		{
			*(mu_avg + i) += *(mu + i);
		}
		for (k = 0; k < CHANNEL; k++)
		{
			for (i = 0; i < sigma_length; i++)
			{
				*(*(sigma_avg + k) + i) += *(*(sigma + k) + i);
			}
			for (i = 0; i < v_length; i++)
			{
				*(*(v_avg + k) + i) += *(*(v + k) + i);
			}
		}
	}
	
			
	for (i = 0; i < CHANNEL; i++)
	{
		*(pi + i) = *(pi_avg + i) / keep_it;
	}
	for (i = 0; i < bigx_cols * CHANNEL; i++)
	{
		*(mu + i) = *(mu_avg + i) / keep_it;
	}
	for (k = 0; k < CHANNEL; k++)
	{
		for (i = 0; i < sigma_length; i++)
		{
			*(*(sigma + k) + i) = *(*(sigma_avg + k) + i) / keep_it;
		}
		for (i = 0; i < v_length; i++)
		{
			*(*(v + k) + i) = *(*(v_avg + k) + i) / keep_it;
		}
	}	
	
	getdenom(intensity, zdenom, bigx, mu, pi, v, sigma, nreads, ncycles, bigx_rows, bigx_cols);

	
	/*
	Destroy all
	*/
	
	for (i = 0; i < CHANNEL * nreads; i++)
	{
		free(*(covgamma + i));
		free(*(index_group + i));
	}
	for (i = 0; i < CHANNEL; i++)
	{
		free(*(sigma_avg + i));
		free(*(v_avg + i));
		free(*(egamma + i));
	}
	free(egamma);
	free(covgamma);
	free(index_group);
	free(sigma_avg);
	free(v_avg);
	free(mu_avg);
	free(pi_avg);
	free(kfn);
	free(bigx);
	free(delta);
	return 0;
}

int init_index_group(int** index_group, int* delta, size_t nreads, size_t ncycles)
{
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t temp = 0;
	
	while (i < nreads)
	{
		while (k < CHANNEL)
		{
			temp = 0;
			/*The first element stores the length of each array*/
			*(index_group + i * CHANNEL + k) = malloc(sizeof(int) * (ncycles * CHANNEL + 1));
			
			while (j < ncycles)
			{
				if (*(delta + i * CHANNEL * ncycles + j * CHANNEL + k) == 1)
				{
					
					*(*(index_group + i * CHANNEL + k) + temp + 1) = j * 4;
					*(*(index_group + i * CHANNEL + k) + temp + 2) = j * 4 + 1;
					*(*(index_group + i * CHANNEL + k) + temp + 3) = j * 4 + 2;
					*(*(index_group + i * CHANNEL + k) + temp + 4) = j * 4 + 3;
					temp += CHANNEL;
				}
				j++;
			}
			*(*(index_group + i * CHANNEL + k) + 0) = temp;
			j = 0;
			k++;
		}
		k = 0;
		i++;
	}
	i = 0;
	
	return 0;
}

int posterior_means(int** index_group, double** egamma, double** covgamma, double** v, double** sigma, double* intensity, double* bigx, double* mu, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols)
{
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t bigv_length = bigx_rows * bigx_rows;
	size_t bigv_scale = bigx_rows;
	double* bigv = malloc(sizeof(double) * bigv_length);
	double* big_sigma12 = malloc(sizeof(double) * bigx_rows * bigx_cols);;
	double* big_sigma22 = malloc(sizeof(double) * bigv_length);
	double* sigma12 = 0;
	double* s12s22 = 0;
	double* mu_k = malloc(sizeof(double) * bigx_cols);
	double* bigx_part = malloc(sizeof(double) * bigx_rows);
	
	
	for (i = 0; i < CHANNEL; i++)
	{
		for (j = 0; j < bigv_length; j++)
		{
			*(bigv + j) = 0.0;
		}
		for (j = 0; j < ncycles; j++)
		{
			size_t temp = j * bigv_scale * CHANNEL + j * CHANNEL;
			for (k = 0; k < CHANNEL; k++)
			{
				*(bigv + temp + k) = *(*(v + i) + k);
				*(bigv + temp + bigv_scale + k) = *(*(v + i) + CHANNEL + k);
				*(bigv + temp + bigv_scale * 2 + k) = *(*(v + i) + CHANNEL * 2 + k);
				*(bigv + temp + bigv_scale * 3 + k) = *(*(v + i) + CHANNEL * 3 + k);
			}
		}
		for (j = 0; j < bigv_length; j++)
		{
			*(big_sigma22 + j) = *(bigv + j);
		}
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bigx_cols, bigx_rows, bigx_cols, 1, *(sigma + i), bigx_cols, bigx, bigx_cols, 0, big_sigma12, bigx_rows);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bigx_rows, bigx_rows, bigx_cols, 1, bigx, bigx_cols, big_sigma12, bigx_rows, 1, big_sigma22, bigx_rows);
		
		for (j = 0; j < bigx_cols; j++)
		{
			*(mu_k + j) = *(mu + j * CHANNEL + i);
		}
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bigx_rows, 1, bigx_cols, 1, bigx, bigx_cols, mu_k, 1, 0, bigx_part, 1);
		
		for (j = 0; j < nreads; j++)
		{
			size_t temp = *(*(index_group + j * CHANNEL + i) + 0);
			for (k = 0; k < bigx_cols; k++)
			{
				*(*(egamma + i) + k * nreads + j) = *(mu_k + k);
			}
			for (k = 0; k < bigx_cols * bigx_cols; k++)
			{
				*(*(covgamma + j * CHANNEL + i) + k) = *(*(sigma + i) + k);
			}
			if (temp > 0)
			{
				int* ipiv = malloc(sizeof(int) * temp);
				double* big_sigma22_part = malloc(sizeof(double) * temp * temp);
				double* egamma_part = malloc(sizeof(double) * temp);
				double* mu_k_temp = malloc(sizeof(double) * bigx_cols);
				
				
				for (k = 0; k < bigx_cols; k++)
				{
					*(mu_k_temp + k) = *(mu_k + k);
				}
		
				for (k = 0; k < temp; k++)
				{
					*(egamma_part + k) = *(intensity + j * bigx_rows + *(*(index_group + j * CHANNEL + i) + k + 1)) - *(bigx_part + *(*(index_group + j * CHANNEL + i) + k + 1));
				}
				
				sigma12 = malloc(sizeof(double) * bigx_cols * temp);
				s12s22 = malloc(sizeof(double) * bigx_cols * temp);
				fetch_cols(big_sigma12, bigx_cols, bigx_rows, *(index_group + j * CHANNEL + i), sigma12);
				fetch_elems(big_sigma22, bigx_rows, bigx_rows, *(index_group + j * CHANNEL + i), *(index_group + j * CHANNEL + i), big_sigma22_part);
				LAPACKE_dgetrf(CblasRowMajor, temp, temp, big_sigma22_part, temp, ipiv);
				LAPACKE_dgetri(CblasRowMajor, temp, big_sigma22_part, temp, ipiv);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bigx_cols, temp, temp, 1, sigma12, temp, big_sigma22_part, temp, 0, s12s22, temp);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bigx_cols, bigx_cols, temp, -1, s12s22, temp, sigma12, temp, 1, *(covgamma + j * CHANNEL + i), bigx_cols);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bigx_cols, 1, temp, 1, s12s22, temp, egamma_part, 1, 1, mu_k_temp, 1);
				
				for (k = 0; k < bigx_cols; k++)
				{
					*(*(egamma + i) + k * nreads + j) = *(mu_k_temp + k);
				}
				
				free(big_sigma22_part);
				free(egamma_part);
				free(mu_k_temp);
				free(sigma12);
				free(s12s22);
			}
			
		}

	}

	free(mu_k);
	free(bigx_part);
	free(bigv);
	free(big_sigma12);
	free(big_sigma22);
	
	return 0;
}

int bigx_row_major(double* bigx, size_t ncycles, size_t nknots, double* knots)
{
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	size_t cols = CHANNEL * (2 + nknots);
	for (i = 0; i < ncycles; i++)
	{
		for (j = 0; j < CHANNEL; j++)
		{
			*(bigx + i * CHANNEL * cols + j * (cols + 1)) = 1;
			*(bigx + i * CHANNEL * cols + j * (cols + 1) + CHANNEL) = i;
			for (k = 0; k < nknots; k++)
			{
				*(bigx + i * CHANNEL * cols + j * (cols + 1) + (k + 2) * CHANNEL) = i > *(knots + k) ? i - *(knots + k) : 0;
			}
		}	
	}
	return 0;
}

int update_pi(int* delta, double* pi, size_t nreads, size_t ncycles)
{
	size_t temp = nreads * ncycles;
	int sums = 0;
	size_t i = 0;
	size_t j = 0;
	for (i = 0; i < CHANNEL; i++)
	{
		sums = 0;
		for (j = 0; j < temp; j++)
		{
			sums += *(delta + j * CHANNEL + i);
		}
		*(pi + i) = ((double)sums)/((double)(temp));
	}
	return 0;
}

int update_mu(double** egamma, double* mu, size_t nreads, size_t bigx_cols)
{
	double* temp_mu = malloc(sizeof(double) * bigx_cols * CHANNEL);
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	for (k = 0; k < CHANNEL; k++)
	{
		for (i = 0; i < bigx_cols; i++)
		{
			*(temp_mu + i * CHANNEL + k) = 0;
			for (j = 0; j < nreads; j++)
			{
				*(temp_mu + i * CHANNEL + k) += *(*(egamma + k) + i * nreads + j);
			}
			*(mu + i * CHANNEL + k) = *(temp_mu + i * CHANNEL + k) / ((double)nreads);
		}
	}
	free(temp_mu);
	return 0;
}

int update_sigma(double** sigma, double** covgamma, double** egamma, double* mu, size_t nreads, size_t bigx_cols)
{
	size_t t = 0;
	size_t i = 0;
	size_t j = 0;
	size_t k = 0;
	for (k = 0; k < CHANNEL; k++)
	{
		double* temp = malloc(sizeof(double) * nreads * bigx_cols);
		double coef = ((double)1)/((double)nreads);
		for (i = 0; i < bigx_cols; i++)
		{
			for (j = 0; j < nreads; j++)
			{
				*(temp + i * nreads + j) = *(*(egamma + k) + i * nreads + j) - *(mu + i * CHANNEL + k);
			}
		}
		for (i = 0; i < bigx_cols; i++)
		{
			for (j = 0; j < bigx_cols; j++)
			{
				*(*(sigma + k) + i * bigx_cols + j) = 0;
				for (t = 0; t < nreads; t++)
				{
					*(*(sigma + k) + i * bigx_cols + j) += *(*(covgamma + t * CHANNEL + k) + i * bigx_cols + j);
				}
			}
		}
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bigx_cols, bigx_cols, nreads, coef, temp, nreads, temp, nreads, coef, *(sigma + k), bigx_cols);
		free(temp);
	}
	return 0;
}

int update_v(double** egamma, double** covgamma, double** v, double* intensity, double* bigx, int* delta, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols)
{
	int value = 1;
	int col_index = 0;
	int delta_sum = 0;
	size_t k = 0;
	size_t j = 0;
	size_t i = 0;
	size_t t = 0;
	double* etemp = 0;
	double* egamma_part = 0;
	double* covtemp = malloc(sizeof(double) * bigx_cols * bigx_cols);
	double* vtemp = malloc(sizeof(double) * CHANNEL * bigx_cols);
	int* selects = malloc(sizeof(int) * (nreads + 1));

	for (k = 0; k < CHANNEL; k++)
	{
		printf("k=%d\n",k);
		for (i = 0; i < CHANNEL * CHANNEL; i++)
		{
			*(*(v + k) + i) = 0;
		}
		for (j = 0; j < ncycles; j++)
		{
			printf("j=%d\n",j);
			col_index = j * CHANNEL + k;
			select_delta(delta, selects, nreads, ncycles, col_index, value);
			printf("%d\n", *(selects + 0));
			if (*(selects + 0) > 0)
			{
				egamma_part = malloc(sizeof(double) * bigx_cols * *(selects + 0));
			etemp = malloc(sizeof(double) * CHANNEL * *(selects + 0));
			for (i = 0; i < bigx_cols * bigx_cols; i++)
			{
				*(covtemp + i) = 0;
			}
			for (i = 0; i < CHANNEL; i++)
			{
				for (t = 0; t < *(selects + 0); t++)
				{
					*(etemp + i * *(selects + 0) + t) = *(intensity + *(selects + t + 1) * bigx_rows + j * CHANNEL + i);
				}
			}
			fetch_cols(*(egamma + k), bigx_cols, nreads, selects, egamma_part);
			printf("0\n");
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, CHANNEL, *(selects + 0), bigx_cols, -1, bigx + j * CHANNEL * bigx_cols, bigx_cols, egamma_part, *(selects + 0), 1, etemp, *(selects + 0));
			for (i = 0; i < bigx_cols * bigx_cols; i++)
			{
				for (t =  0; t < *(selects + 0); t++)
				{
					*(covtemp + i) += *(*(covgamma + *(selects + t + 1) * CHANNEL + k) + i);					
				}
			}
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, CHANNEL, bigx_cols, bigx_cols, 1, bigx + j * CHANNEL * bigx_cols, bigx_cols, covtemp, bigx_cols, 0, vtemp, bigx_cols);
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, CHANNEL, CHANNEL, bigx_cols, 1, vtemp, bigx_cols, bigx + j * CHANNEL * bigx_cols, bigx_cols, 1, *(v + k), CHANNEL);
			print_matrix_rowmajor("etemp", CHANNEL, *(selects + 0), etemp, *(selects + 0));
			print_matrix_rowmajor("vk", CHANNEL, CHANNEL, *(v+k), CHANNEL);
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, CHANNEL, CHANNEL, *(selects + 0), 1, etemp, *(selects + 0), etemp, *(selects + 0), 1, *(v + k), CHANNEL);
			printf("4\n");
			free(etemp);
			free(egamma_part);
			}
		}
		delta_sum = 0;
		for (i = 0; i < nreads; i++)
		{
			for (j = 0; j < ncycles; i++)
			{
				delta_sum += *(delta + i * bigx_rows + j * CHANNEL + k);
			}
		}
		for (i = 0; i < CHANNEL * CHANNEL; i++)
		{
			*(*(v + k) + i) /= delta_sum;
		}
	}

	free(vtemp);
	free(covtemp);
	free(selects);
	return 0;
}

int getdenom(double* intensity, double* denom, double* bigx, double* mu, double* pi, double** v, double** sigma, size_t nreads, size_t ncycles, size_t bigx_rows, size_t bigx_cols)
{
	size_t k = 0;
	size_t j = 0;
	size_t i = 0;
	size_t t = 0;
	size_t intensity_length = nreads * ncycles;
	double d = 0.0;
	double* mu_k = malloc(sizeof(double) * bigx_cols);
	double* mean = malloc(sizeof(double) * CHANNEL);
	double* tempv = malloc(sizeof(double) * CHANNEL * CHANNEL);
	double* temp_sigma = malloc(sizeof(double) * CHANNEL * bigx_cols);
	double* int_temp = malloc(sizeof(double) * CHANNEL);
	double* final_temp = malloc(sizeof(double) * CHANNEL);
	double* final = malloc(sizeof(double) * 1);
	int* ipiv = malloc(sizeof(int) * CHANNEL);
	for(k = 0; k < CHANNEL; k++)
	{
		for (i = 0; i < bigx_cols; i++)
		{
			*(mu_k + i) = *(mu + i * CHANNEL + k);
		}
		
		for (j = 0; j < ncycles; j++)
		{
			for (i = 0; i < CHANNEL * CHANNEL; i++)
			{
				*(tempv + i) = *(*(v + k) + i);
			}
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, CHANNEL, 1, bigx_cols, 1, bigx + j * CHANNEL * bigx_cols, bigx_cols, mu_k, 1, 0, mean, 1);
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, CHANNEL, bigx_cols, bigx_cols, 1, bigx + j * CHANNEL * bigx_cols, bigx_cols, *(sigma + k), bigx_cols, 0, temp_sigma, bigx_cols);
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, CHANNEL, CHANNEL, bigx_cols, 1, temp_sigma, bigx_cols, bigx + j * CHANNEL * bigx_cols, bigx_cols, 1, tempv, CHANNEL);
			d = pow(det(tempv, CHANNEL, ipiv), -0.5) * pow(2 * PI, -CHANNEL / 2);
			LAPACKE_dgetri(CblasRowMajor, CHANNEL, tempv, CHANNEL, ipiv);
			for (i = 0; i < nreads; i++)
			{
				for (t = 0; t < CHANNEL; t++)
				{
					*(int_temp + t) = *(intensity + i * bigx_rows + j * CHANNEL + t) - *(mean + t);
				}
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 1, CHANNEL, CHANNEL, 1, int_temp, 1, tempv, CHANNEL, 0, final_temp, CHANNEL);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 1, CHANNEL, 1, final_temp, CHANNEL, int_temp, 1, 0, final, 1);
				*(denom + i * bigx_rows + j * ncycles + k) = d * exp(*(final + 0) * (-0.5));
			}
		}
	}
	
	for (i = 0; i < intensity_length; i++)
	{
		double temp_sum = SMALL * CHANNEL;
		for (k = 0; k < CHANNEL; k++)
		{
			*(denom + k + i * CHANNEL) += SMALL;
		}
		for (k = 0; k < CHANNEL; k++)
		{
			*(denom + k + i * CHANNEL) *= *(pi + k);
		}
		for (k = 0; k < CHANNEL; k++)
		{
			temp_sum += *(denom + k + i * CHANNEL);
		}
		for (k = 0; k < CHANNEL; k++)
		{
			*(denom + k + i * CHANNEL) /= temp_sum;
		}
	}
	
	free(ipiv);
	free(final);
	free(final_temp);
	free(int_temp);
	free(temp_sigma);
	free(tempv);
	free(mean);
	free(mu_k);
	return 0;
}

int update_delta(double* denom, int* delta, double q, size_t nreads, size_t ncycles)
{
	size_t i = 0;
	size_t k = 0;
	size_t length =  nreads * ncycles;
	for (i = 0; i < length; i++)
	{
		for (k = 0; k < CHANNEL; k++)
		{
			*(delta + i * CHANNEL + k) = 0;
		}
		*(delta + rmultinom_array(denom, i * CHANNEL, (i + 1) * CHANNEL)) = 1;
	}
	return 0;
}

double det(double* matrix, size_t length, int* ipiv)
{
	double result = 1;
	size_t i = 0;
	int sum = 0;
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, length, length, matrix, length, ipiv);
	for (i = 0; i < length; i++)
	{
		result *= *(matrix + i * length + i);
		if ((i + 1) != *(ipiv + i))
		{
			sum++;
		}
	}
	if (sum % 2 != 0)
	{
		result *= -1;
	} 
	return result;
}