#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "cephes.h"
#include "rats.h"

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
I T E M ---- L I N E A R  C O M P L E X I T Y  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//data: input binary sequence
//offset: start (may be 0)
//bits: sequence length

int rats_Berlekamp_Massey(unsigned char *data, int offset, int bits)
{
	int i, j, N, n, L, m, d, U, b, r;

#define Berlekamp_Massey_Max_Sequence_Unit		32
	unsigned int s[Berlekamp_Massey_Max_Sequence_Unit];
	unsigned int B[Berlekamp_Massey_Max_Sequence_Unit];
	unsigned int C[Berlekamp_Massey_Max_Sequence_Unit];
	unsigned int P[Berlekamp_Massey_Max_Sequence_Unit];

	if (bits > 8 * sizeof(unsigned int)*Berlekamp_Massey_Max_Sequence_Unit)
		return -1;
#undef	Berlekamp_Massey_Max_Sequence_Unit

	n = bits;
	U = (n + 31) / 32;

	//Initialization.
	L = 0;	m = -1;	d = 0;
	for (i = 0; i < U; i++)
	{
		s[i] = 0;
		B[i] = 0;
		C[i] = 0;
		P[i] = 0;
	}
	B[0] = 1;
	C[0] = 1;

	//Make sequence s[n-1],s[n-2],...,s[1],s[0]
	for (i = 0; i < n; i++)
	{
		s[i / 32] = (s[i / 32] << 1) + GET_EPSILON(data, offset + n - i -1);
	}

	/* DETERMINE LINEAR COMPLEXITY */
	for (N = 0; N < n; N++)
	{
		d = GET_EPSILON(data, offset + N);//C.bit0 == 1 forever
		for (i = 1; i <= L; i++)
			d += ((C[i / 32] >> (31 - i)) & 1) * GET_EPSILON(data, offset + N - i);
		if (d = d % 2)
		{
			//P(D) = B(D) · D^(N−m).
			//We know N - m >= 1
			r = (N - m) % 8;
			if (b = (N - m) / 8)
			{
				for (j = U - 1; j > b; j--)
					P[j] = (B[j - b] >> r) | (B[j - b - 1] << (32 - r));
				P[j--] = B[0] >> r;
				while (j >= 0)	P[j--] = 0;
			}
			else
			{
				for (j = U - 1; j > 0; j--)
					P[j] = (B[j] >> r) | (B[j - 1] << (32 - r));
				P[0] = B[0] >> r;
			}
			if (L <= N / 2)
			{
				L = N + 1 - L;
				m = N;
				for (i = 0; i < U; i++)
				{
					B[i] = C[i];
					C[i] ^= P[i];
				}
			}
			else
			{
				for (i = 0; i < U; i++)
				{
					C[i] ^= P[i];
				}
			}
		}
	}

	return L;
}


int LinearComplexity0(double alpha, unsigned char *data, int bits, int M, LinearComplexity_V *value)
{
	int       i, ii, j, d, N, L, m, N_, sign, K = 6;
	double    p_value, v_obs, Ti, u;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	int       v[7] = { 0, 0, 0, 0, 0, 0, 0 };
	unsigned char  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	time_t s,e;
	float dur = 0;


	N = bits / M;
	if (((B_ = (unsigned char  *)calloc(M, sizeof(unsigned char))) == NULL) ||
		((C = (unsigned char  *)calloc(M, sizeof(unsigned char))) == NULL) ||
		((P = (unsigned char  *)calloc(M, sizeof(unsigned char))) == NULL) ||
		((T = (unsigned char  *)calloc(M, sizeof(unsigned char))) == NULL)) 
	{
		if (B_ != NULL)		free(B_);
		if (C != NULL)		free(C);
		if (P != NULL)		free(P);
		if (T != NULL)		free(T);
		return -1;		//ERROR: Insufficient Memory for Work Space.
	}

	//calculate the theoretical mean
	sign = ((M + 1) & 1) ? -1 : 1;
	printf("sign:%d",sign);
	u = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M) * (M / 3.0 + 2.0 / 9.0);
	sign = (M & 1) ? -1 : 1;
	printf("M:%d , u:%f , sign:%d\n",M,u,sign);


	for (ii = 0; ii<N; ii++)
	{
		s = clock();
		for (i = 0; i<M; i++) {
			B_[i] = 0;
			C[i] = 0;
			T[i] = 0;
			P[i] = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0] = 1;
		B_[0] = 1;

		/* DETERMINE LINEAR COMPLEXITY */
		for (N_ = 0; N_ < M; N_++)
		{
			d = GET_EPSILON(data, ii*M + N_);
			for (i = 1; i <= L; i++)
				d += C[i] * GET_EPSILON(data, ii*M + N_ - i);
			if (d = d % 2)
			{
				for (i = 0; i < M; i++)
				{
					T[i] = C[i];
					P[i] = 0;
				}
				for (j = 0; j < M; j++)
				{
					if (B_[j] == 1)
						P[j + N_ - m] = 1;
				}
				for (i = 0; i < M; i++)
					C[i] = (C[i] + P[i]) % 2;
				if (L <= N_ / 2)
				{
					L = N_ + 1 - L;
					m = N_;
					for (i = 0; i<M; i++)
						B_[i] = T[i];
				}
			}
		}

		//if (L != rats_Berlekamp_Massey(data, ii*M, M))
		//{
		//	L = rats_Berlekamp_Massey(data, ii*M, M);
		//}

		//Calculate Ti, and record the Ti values in v0,…, v6
		Ti = sign * (L - u) + 2.0 / 9.0;

		if (Ti <= -2.5)
			v[0]++;
		else if (Ti <= -1.5)
			v[1]++;
		else if (Ti <= -0.5)
			v[2]++;
		else if (Ti <= 0.5)
			v[3]++;
		else if (Ti <= 1.5)
			v[4]++;
		else if (Ti <= 2.5)
			v[5]++;
		else
			v[6]++;

	e = clock();
	dur = (float)(e-s)/CLOCKS_PER_SEC;
	//printf("%d liner time:%f\n",ii,dur);
	}



	free(B_);
	free(P);
	free(C);
	free(T);

	//Compute P-value
	v_obs = 0.00;
	for (i = 0; i< 7; i++)
	{
		printf("vi: %d\n",v[i]);
		v_obs += pow((double)v[i] - N*pi[i], 2) / (N*pi[i]);
	}
	p_value = cephes_igamc(3.0, v_obs / 2.0);

	printf("obs: %f\n",v_obs);

	//TODO: return value

	if (p_value < alpha)
		return 0;
	return 1;
}

int LinearComplexity1(double alpha, unsigned char *data, int bits, int M, LinearComplexity_V *value){

	int   sign, K = 6;
	double u;
	sign = ((M + 1) & 1) ? -1 : 1;
	u = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M) * (M / 3.0 + 2.0 / 9.0);
	sign = (M & 1) ? -1 : 1;
//	printf("M:%d , u:%f , sign:%d\n",M,u,sign);


	double v_obs = LinearComplexity(alpha,data,bits,M,u,sign,value);
//	printf("v_obs: %f\n",v_obs);
	double p_value = cephes_igamc(3.0, v_obs / 2.0);



	//TODO: return value

	if (p_value < alpha)
		return 0;
	return 1;
}