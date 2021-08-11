#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "cephes.h"
#include "rats.h"



#define MAXNUMOFTEMPLATES				148		/* APERIODIC TEMPLATES: 148=>temp_length=9 */
#define MIN(x,y)   ((x) >  (y)  ? (y)  : (x))

typedef struct para{
	int M;
	int N;
	int m;
	double lambda;
	double varWj;
} gpu_param;

__device__ int GET_EPSILON1(unsigned char*data,int offset){

	return ((data[offset/8]>>(7-(offset)%8))&1);
}

__device__ int MINVALUE(int a,int b){

	return a > b ? b : a;
}

 
__global__ void ReductionSum(double *CHI2_gpu, gpu_param *M_gpu, unsigned char* data, unsigned char*sequence){

	int  i, jj , k , j , match, K = 5;;
	unsigned int  W_obs, nu[6];
	double sum,chi2;
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};

	int N = M_gpu->N;
	int M = M_gpu->M;
	int m = M_gpu->m;
	double lambda = M_gpu->lambda;
	double varWj = M_gpu->varWj;


	unsigned int *Wj = NULL;
	if ( (Wj = (unsigned int*)malloc(N*sizeof(unsigned int))) == NULL ) {
			return ;
		}
	//确定索引
	int threadId = blockIdx.x *blockDim.x + threadIdx.x;  


	if(threadId<MINVALUE(MAXNUMOFTEMPLATES, numOfTemplates[m]))
	{

		//printf("N:%d , M:%d , m:%d , lam:%f , var:%f\n",N,M,m,lambda,varWj);
		sum = 0;
		for ( k=0; k<=K; k++ )
			nu[k] = 0;
		for ( i=0; i<N; i++ ) {
			W_obs = 0;
			for ( j=0; j<M-m+1; j++ ) {
				match = 1;
				for ( k=0; k<m; k++ ) {
					if ( (int)sequence[threadId*m+k] != (int)GET_EPSILON1(data,i*M+j+k) ) {
						match = 0;
						break;
					}
				}
				if ( match == 1 ) {
					W_obs++;
					j += m-1;
				}
			}
			Wj[i] = W_obs;
		}
		sum = 0;
		chi2 = 0.0;                                   /* Compute Chi Square */
		for ( i=0; i<N; i++ ) {
		//	if ( m == 10 )
			//	fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
		//	else
			//	fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
			chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
		}
		CHI2_gpu[threadId] = chi2;
	}

	


}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int
NonOverlappingTemplateMatchingsCUDA(double alpha, unsigned char *data, int bits, int m)
{
    int n = bits;
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, nu[6], *Wj = NULL; 
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	unsigned char		*sequence = NULL;
	unsigned char		sequence1[(MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]))][m] ;
	int i_,k_;


	time_t s,e;

	N = 8;
	M = n/N;

	if ( (Wj = (unsigned int*)malloc(N*sizeof(unsigned int))) == NULL ) {
		return 0;
	}
	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	sprintf(directory, "../templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL|| (((sequence = (unsigned char*)malloc(MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*m*sizeof(unsigned char))) == NULL )))) {
		if ( sequence != NULL )
			free(sequence); 
	}
	else {

		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

	/////////////////////////////////////////////
	for (i_ = 0; i_ < (MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])); i_++){
			for (k_=0; k_<m; k_++ ) {
				fscanf(fp, "%d", &bit);
				sequence1[i_][k_] = bit;
				sequence[i_*m+k_] = bit;
			}
	}
	/////////////////////////////////////////////
	//分配host变量
	gpu_param M_cpu;
	M_cpu.M = M;
	M_cpu.N = N;
	M_cpu.m = m;
	M_cpu.lambda = lambda;
	M_cpu.varWj = varWj;
	double CHI2[MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])] = {0.0};


	//分配device变量
	gpu_param *M_gpu;
	unsigned char *data_gpu;
	double *CHI2_gpu;
	cudaMalloc((void**)&CHI2_gpu,(MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*sizeof(double)));

	unsigned char *sequenceGPU ;
	cudaMalloc((void**)&sequenceGPU,MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*m*sizeof(double));

	cudaMalloc((void**)&M_gpu,sizeof(gpu_param));
	cudaMalloc((void**)&data_gpu,(bits / 8)*sizeof(unsigned char));

	cudaMemcpy(M_gpu, &M_cpu, sizeof(gpu_param), cudaMemcpyHostToDevice);
	cudaMemcpy(data_gpu, data, (bits / 8)*sizeof(unsigned char), cudaMemcpyHostToDevice);
	cudaMemcpy(sequenceGPU, sequence, MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*m*sizeof(double), cudaMemcpyHostToDevice);

	int blocksPerGrid = (MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])+512-1)/512;

	ReductionSum<<<blocksPerGrid,512>>>(CHI2_gpu,M_gpu,data_gpu,sequenceGPU);


	//将结果传回到主机端
	cudaMemcpy(CHI2, CHI2_gpu, MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*sizeof(double), cudaMemcpyDeviceToHost);


		for( jj=0; jj<MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++ ) {
			//printf("CHI2:%f    ",CHI2[jj]);
			p_value = cephes_igamc(N/2.0, CHI2[jj]/2.0);
		
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				printf("\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

		//	printf("%9.6f %f %s %3d\n", CHI2[jj], p_value, p_value < alpha ? "FAILURE" : "SUCCESS", jj);
			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
		}
	}

	if ( sequence != NULL )
		free(sequence);

	free(Wj);

    if ( fp != NULL )
        fclose(fp);
   // printf("p_value: %f\n",p_value);
    if (p_value < alpha)
		return 0;
	return 1;
}
