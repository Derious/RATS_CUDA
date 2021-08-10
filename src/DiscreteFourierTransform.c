#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"

#ifndef SQRT2
#define SQRT2		1.41421356237309504880
#endif // !SQRT2

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void  __ogg_fdrffti(int n, double *wsave, int *ifac);
void  __ogg_fdrfftf(int n, double *X, double *wsave, int *ifac);



int DiscreteFourierTransform(double alpha, unsigned char *data, int bits, DFT_V *value)
{
	double	p_value, T, percentile, N_l, N_0, v_obs;
	double *m = NULL, *X = NULL, *wsave = NULL;
	int		i, count, ifac[15];

	if (((X = (double*)calloc(bits, sizeof(double))) == NULL) ||
		((wsave = (double *)calloc(2 * bits, sizeof(double))) == NULL) ||
		((m = (double*)calloc(bits / 2 + 1, sizeof(double))) == NULL)) 
	{
		if (X != NULL)		free(X);
		if (wsave != NULL)	free(wsave);
		if (m != NULL)		free(m);
		return -1;		//ERROR: Unable to allocate working arrays for the DFT.
	}
	for (i = 0; i<bits; i++)
		X[i] = 2 * GET_EPSILON(data, i) - 1;

	__ogg_fdrffti(bits, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(bits, X, wsave, ifac);	/* APPLY FORWARD FFT */

	m[0] = sqrt(X[0] * X[0]);	    /* COMPUTE MAGNITUDE */
	for (i = 0; i<bits / 2; i++)
		m[i + 1] = sqrt(pow(X[2 * i + 1], 2) + pow(X[2 * i + 2], 2));

	T = sqrt(2.995732274*bits);
	count = 0;				       /* CONFIDENCE INTERVAL */
	for (i = 0; i < bits / 2; i++)
	{
		if (m[i] < T)
			count++;
	}
	percentile = (double)count / (bits / 2) * 100;
	N_l = (double)count;       /* number of peaks less than h = sqrt(3*bits) */
	N_0 = (double) 0.95*bits / 2.0;
	v_obs = (N_l - N_0) / sqrt(bits / 4.0 * 0.95 * 0.05);
	p_value = erfc(fabs(v_obs) / SQRT2);

	free(X);
	free(wsave);
	free(m);

	//TODO: return value

	if (p_value < alpha)
		return 0;
	return 1;
}
