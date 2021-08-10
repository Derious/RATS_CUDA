#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"

#ifndef SQRT2
#define SQRT2		1.41421356237309504880
#endif // !SQRT2

#ifndef GET_UINT24_BE
#define GET_UINT24_BE(d, n)	(((unsigned int)(d)[n] << 16) | ((unsigned int)(d)[(n)+1] << 8) | ((unsigned int)(d)[(n)+2]))
#endif // !GET_UINT24_BE

#ifndef GET_2EX_SUB_ONE
#define GET_2EX_SUB_ONE(x)	(((unsigned int)1 << (x)) - 1)
#endif // !GET_2EX_SUB_ONE

#ifndef GET_SEGMENT_U24
//Get 'L' bits segment of byte string 'd'(from 's' to 's + L - 1', s >= 0)
#define GET_SEGMENT_U24(d, s, L)	( (GET_UINT24_BE(d, (s) / 8) >> (24 - L - (s) % 8))	& GET_2EX_SUB_ONE(L) )
#endif // !GET_SEGMENT_U24

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                 M A U R E R - U N I V E R S A L - T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


int MaurerUniversal(double alpha, unsigned char *data, int bits, MaurerUniversal_V *value)
{
	int		i, j, p, L, Q, K;
	double	p_value, v_obs, sigma, fn, sum, c;
	int	*T;
	double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
		8.1764248, 9.1723243, 10.170032, 11.168765,
		12.168070, 13.167693, 14.167488, 15.167379 };
	double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
		3.401, 3.410, 3.416, 3.419, 3.421 };

	/* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
	 * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
	 * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if (bits <  387840)				L = 5;
	else if (bits <  904960)		L = 6;
	else if (bits <  2068480)		L = 7;
	else if (bits <  4654080)		L = 8;
	else if (bits <  10342400)		L = 9;
	else if (bits <  22753280)		L = 10;
	else if (bits <  49643520)		L = 11;
	else if (bits <  107560960)		L = 12;
	else if (bits <  231669760)		L = 13;
	else if (bits <  496435200)		L = 14;
	else if (bits <  1059061760)	L = 15;
	else							L = 16;

	p = (int)pow(2, L);
	Q = 10 * p;
	K = bits / L - Q;	 		    /* BLOCKS TO TEST */

	if ((L < 6) || (L > 16))
		return -1;				//ERROR:  bits IS OUT OF RANGE.
	if ((T = (int *)calloc(p, sizeof(int))) == NULL)
		return -2;				//ERROR: Unable to allocate T.

	/* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
	c = 0.7 - 0.8 / (double)L + (4 + 32 / (double)L)*pow(K, -3 / (double)L) / 15;
	sigma = c * sqrt(variance[L] / (double)K);
	
	/* INITIALIZE TABLE */
	memset(T, 0, p * sizeof(int));
	for (i = 1; i <= Q; i++) {
		j = GET_SEGMENT_U24(data, (i - 1) * L, L);
		T[j] = i;
	}
	/* PROCESS BLOCKS */
	sum = 0.0;
	for (i = Q + 1; i <= Q + K; i++) {
		j = GET_SEGMENT_U24(data, (i - 1) * L, L);
		sum += log(i - T[j]) / log(2);
		T[j] = i;
	}
	fn = (double)(sum/(double)K);

	v_obs = fabs(fn - expected_value[L]) / sigma;
	p_value = erfc(v_obs / SQRT2);

	free(T);

	//TODO: return value

	if (p_value < alpha)
		return 0;
	return 1;
}