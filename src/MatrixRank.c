#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cephes.h"
#include "matrix.h"
#include "rats.h"

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

int MatrixRank(double alpha, unsigned char *data, int bits, int m, MatrixRank_V *value)
{
	int			N, i,j, k;
	double		p_value, v_obs, arg1, p_32, p_31, p_30, R;
	unsigned char	**matrix;
	int			F_32, F_31, F_30;
	int			M = 32, Q = 32;

	N = bits / (32 * 32);
	if (N == 0)
		return -1;		//ERROR: Insuffucient # Of Bits To Define An 32x32  Matrix

	/* COMPUTE PROBABILITIES */
	//r = 32;
	//product = 1;
	//for (i = 0; i <= r - 1; i++)
	//	product *= ((1.e0 - pow(2, i - 32))*(1.e0 - pow(2, i - 32))) / (1.e0 - pow(2, i - r));
	//p_32 = pow(2, r*(32 + 32 - r) - 32 * 32) * product;

	//r = 31;
	//product = 1;
	//for (i = 0; i <= r - 1; i++)
	//	product *= ((1.e0 - pow(2, i - 32))*(1.e0 - pow(2, i - 32))) / (1.e0 - pow(2, i - r));
	//p_31 = pow(2, r*(32 + 32 - r) - 32 * 32) * product;

	//p_30 = 1 - (p_32 + p_31);

	p_32 = 0.2888;
	p_31 = 0.5776;
	p_30 = 0.1336;

	matrix = create_matrix(M, Q);

	F_32 = 0;
	F_31 = 0;
	for (k = 0; k<N; k++) 
	{/* FOR EACH 32x32 MATRIX   */
		for (i = 0; i < M; i++)
			for (j = 0; j < Q; j++)
				matrix[i][j] = GET_EPSILON(data, k*(M*Q) + i*M + j);

		R = computeRank(32, 32, matrix);
		if (R == 32)	F_32++;			/* DETERMINE FREQUENCIES */
		if (R == 31)	F_31++;
	}
	F_30 = N - F_32 - F_31;

	v_obs = (pow(F_32 - N*p_32, 2) / (double)(N*p_32) +
		pow(F_31 - N*p_31, 2) / (double)(N*p_31) +
		pow(F_30 - N*p_30, 2) / (double)(N*p_30));

	arg1 = -v_obs / 2.0;

	p_value = exp(arg1);

	delete_matrix(M, matrix);

	//TODO: return value & ¿ÉÓÅ»¯

	if (p_value < alpha)
		return 0;
	return 1;
}
