/* got rid of unused 'k' */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cephes.h"  
#include "rats.h"

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

int LongestRun(double alpha, unsigned char *data, int bits, LongstRun_V *value)
{
	double			p_value, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M, V[7];
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

	if (bits < 128) {
		return -1;		//ERROR: bits is too short.
	}
	else if (bits < 6272) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.2148;
		pi[1] = 0.3672;
		pi[2] = 0.2305;
		pi[3] = 0.1875;
	}
	else if (bits < 750000) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174;
		pi[1] = 0.2430;
		pi[2] = 0.2493;
		pi[3] = 0.1752;
		pi[4] = 0.1027;
		pi[5] = 0.1124;
	}
	else {
		K = 6;
		M = 10000;
		V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}

	N = bits / M;
	for (i = 0; i < N; i++) {
		v_n_obs = 0;
		run = 0;
		//We know M%8 == 0
		//TODO: 可优化
		for (j = 0; j < M; j++) {
			if (GET_EPSILON(data, i*M + j) == 1) {
				run++;
				if (run > v_n_obs)
					v_n_obs = run;
			}
			else
				run = 0;
		}

		if (v_n_obs <= V[0])
			nu[0]++;
		else if (v_n_obs >= V[K])
			nu[K]++;
		else
			nu[v_n_obs - V[0]]++;
	}

	chi2 = 0.0;
	for (i = 0; i <= K; i++)
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	p_value = cephes_igamc((double)(K / 2.0), chi2 / 2.0);

	//TODO: return value & 可优化
	if (value)
	{
		value->p_value = p_value;
	}

	if (p_value < alpha)
		return 0;
	return 1;
}
