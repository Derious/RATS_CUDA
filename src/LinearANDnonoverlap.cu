#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rats.h"
#include "cephes.h"

#define BITS	1000000
#define ROUND	1

typedef struct para{
	int M;
	double u;
	int sign;
	int N;
} gpu_param;

#define MAXNUMOFTEMPLATES				148		/* APERIODIC TEMPLATES: 148=>temp_length=9 */
#define MIN(x,y)   ((x) >  (y)  ? (y)  : (x))

typedef struct para1{
	int M;
	int N;
	int m;
	double lambda;
	double varWj;
} non_param;

__device__ int GET_EPSILON_line(unsigned char*data,int offset){

	return ((data[offset/8]>>(7-(offset)%8))&1);
}
__device__ int MINVALUE1(int a,int b){

	return a > b ? b : a;
}

// __device__ int rats_Berlekamp_Massey(unsigned char *data, int offset, int bits)
// {
// 	int i, j, N, n, L, m, d, U, b, r;

// #define Berlekamp_Massey_Max_Sequence_Unit		32
// 	unsigned int s[Berlekamp_Massey_Max_Sequence_Unit];
// 	unsigned int B[Berlekamp_Massey_Max_Sequence_Unit];
// 	unsigned int C[Berlekamp_Massey_Max_Sequence_Unit];
// 	unsigned int P[Berlekamp_Massey_Max_Sequence_Unit];

// 	if (bits > 8 * sizeof(unsigned int)*Berlekamp_Massey_Max_Sequence_Unit)
// 		return -1;
// #undef	Berlekamp_Massey_Max_Sequence_Unit

// 	n = bits;
// 	U = (n + 31) / 32;

// 	//Initialization.
// 	L = 0;	m = -1;	d = 0;
// 	for (i = 0; i < U; i++)
// 	{
// 		s[i] = 0;
// 		B[i] = 0;
// 		C[i] = 0;
// 		P[i] = 0;
// 	}
// 	B[0] = 1;
// 	C[0] = 1;

// 	//Make sequence s[n-1],s[n-2],...,s[1],s[0]
// 	for (i = 0; i < n; i++)
// 	{
// 		s[i / 32] = (s[i / 32] << 1) + GET_EPSILON(data, offset + n - i -1);
// 	}

// 	/* DETERMINE LINEAR COMPLEXITY */
// 	for (N = 0; N < n; N++)
// 	{
// 		d = GET_EPSILON(data, offset + N);//C.bit0 == 1 forever
// 		for (i = 1; i <= L; i++)
// 			d += ((C[i / 32] >> (31 - i)) & 1) * GET_EPSILON(data, offset + N - i);
// 		if (d = d % 2)
// 		{
// 			//P(D) = B(D) · D^(N−m).
// 			//We know N - m >= 1
// 			r = (N - m) % 8;
// 			if (b = (N - m) / 8)
// 			{
// 				for (j = U - 1; j > b; j--)
// 					P[j] = (B[j - b] >> r) | (B[j - b - 1] << (32 - r));
// 				P[j--] = B[0] >> r;
// 				while (j >= 0)	P[j--] = 0;
// 			}
// 			else
// 			{
// 				for (j = U - 1; j > 0; j--)
// 					P[j] = (B[j] >> r) | (B[j - 1] << (32 - r));
// 				P[0] = B[0] >> r;
// 			}
// 			if (L <= N / 2)
// 			{
// 				L = N + 1 - L;
// 				m = N;
// 				for (i = 0; i < U; i++)
// 				{
// 					B[i] = C[i];
// 					C[i] ^= P[i];
// 				}
// 			}
// 			else
// 			{
// 				for (i = 0; i < U; i++)
// 				{
// 					C[i] ^= P[i];
// 				}
// 			}
// 		}
// 	}

// 	return L;
// }


 
__global__ void ReductionSum(int* WJ, double *CHI2_gpu, non_param *Mnon_gpu,unsigned char*sequence,double *d_partial_sum, gpu_param *M_gpu, unsigned char* data)
{
//printf("11111111111111111");
	//申请共享内存，存在于每个block中 
	int M = M_gpu->M;
	double u = M_gpu->u;
	double Ti;
	int sign = M_gpu->sign;
	int N = M_gpu->N;
	//printf("M = %d\n",M);

	//__shared__ int v[7];
	//确定索引
	int threadId = blockIdx.y * blockDim.y + threadIdx.y;  
   // int threadId = blockIdx.x *blockDim.x + threadIdx.x;  
	int blockID = blockIdx.x;
	int tid = threadIdx.x;
//	int i = threadIdx.x + blockIdx.x * blockDim.x;
	//int tid = threadIdx.x;
	int L , m , d , N_;
	int i , j;
//printf("threadid:%d,blcokid:%d,tid:%d\n",threadId,blockID,tid);
if(threadId<N && blockIdx.x == 0 && threadIdx.x == 0){
	//printf("threadID:%d\n",threadId);
	unsigned char *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	T = (unsigned char*)malloc(M*sizeof(unsigned char));
	P = (unsigned char*)malloc(M*sizeof(unsigned char));
	B_ = (unsigned char*)malloc(M*sizeof(unsigned char));
	C = (unsigned char*)malloc(M*sizeof(unsigned char));
	for (int i = 0; i<M; i++) {
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
			d = GET_EPSILON_line(data, threadId*M + N_);
			//printf("d: %d\n",d);
			for ( j = 1; j <= L; j++)
				d += C[j] * GET_EPSILON_line(data, threadId*M + N_ - j);
			if (d = d % 2)
			{
				for ( j = 0; j < M; j++)
				{
					T[j] = C[j];
					P[j] = 0;
				}
				for ( j = 0; j < M; j++)
				{
					if (B_[j] == 1)
						P[j + N_ - m] = 1;
				}
				for ( j = 0; j < M; j++)
					C[j] = (C[j] + P[j]) % 2;
				if (L <= N_ / 2)
				{
					L = N_ + 1 - L;
					m = N_;
					for ( j = 0; j<M; j++)
						B_[j] = T[j];
				}
			}
		}
		//L = rats_Berlekamp_Massey(data, threadId*M, M);
		Ti = sign * (L - u) + 2.0 / 9.0;
		d_partial_sum[threadId] = Ti;
		//printf("threadId:%d Ti: %f\n",threadId,Ti);


	//printf("M:%d u:%f sign:%d Ti:%f\n",M,u,sign,Ti);
	// free(B_);
	// free(P);
	// free(C);
	// free(T);
}

int  jj , k  , match, K = 5;;
unsigned int  W_obs, nu[6];
double sum,chi2;
int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
    2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};

N = Mnon_gpu->N;
M = Mnon_gpu->M;
m = Mnon_gpu->m;
double lambda = Mnon_gpu->lambda;
double varWj = Mnon_gpu->varWj;



if(tid<MINVALUE1(MAXNUMOFTEMPLATES, numOfTemplates[m]) && blockIdx.y == 0 && threadIdx.y == 0)
{

    //printf("N:%d , M:%d , m:%d , lam:%f , var:%f\n",N,M,m,lambda,varWj);
    //sum = 0;
    // for ( k=0; k<=K; k++ )
    // 	nu[k] = 0;
    // for ( i=0; i<N; i++ ) {
        W_obs = 0;
        for ( j=0; j<M-m+1; j++ ) {
            match = 1;
            for ( k=0; k<m; k++ ) {
                if ( (int)sequence[tid*m+k] != (int)GET_EPSILON_line(data,blockID*M+j+k) ) {
                    match = 0;
                    break;
                }
            }
            if ( match == 1 ) {
                W_obs++;
                j += m-1;
            }
        }
        //Wj[tid*N+blockID] = W_obs;
        WJ[tid*N+blockID] = W_obs;
  //  printf("tid:%d bid:%d wj:%d\n",tid,blockID,WJ[tid*N+blockID]);
}
}
 
double LinearANDnonoverlap(double alpha, unsigned char *data, int bits, int M, double u, int sign, int m, LinearComplexity_V *value)
{

    /////////////////////////////////////
    //非重复模板匹配所需要的参数
    /////////////////////////////////////
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
	int				i, j, jj, k, match, SKIP, M_, N, K = 5;
	char			directory[100];
	unsigned char		*sequence = NULL;
	unsigned char		sequence1[(MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]))][m] ;
	int i_,k_;

	N = 8;
	M_ = n/N;

	if ( (Wj = (unsigned int*)malloc(N*sizeof(unsigned int))) == NULL ) {
		return 0;
	}
	lambda = (M_-m+1)/pow(2, m);
	varWj = M_*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
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
		non_param Mnon_cpu;
		Mnon_cpu.M = M_;
		Mnon_cpu.N = N;
		Mnon_cpu.m = m;
		Mnon_cpu.lambda = lambda;
		Mnon_cpu.varWj = varWj;
		double CHI2[MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])] = {0.0};
		int WJ[148*8] = {0};

	//分配device变量

	


    non_param *Mnon_gpu;
    //unsigned char *data_gpu;
    double *CHI2_gpu;
    cudaMalloc((void**)&CHI2_gpu,(MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*sizeof(double)));

    unsigned char *sequenceGPU ;
    cudaMalloc((void**)&sequenceGPU,MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*m*sizeof(double));

    cudaMalloc((void**)&Mnon_gpu,sizeof(non_param));
    //cudaMalloc((void**)&data_gpu,(bits / 8)*sizeof(unsigned char));

    int * WJ_gpu;
    cudaMalloc((void**)&WJ_gpu,148*8*sizeof(int));

    cudaMemcpy(Mnon_gpu, &Mnon_cpu, sizeof(non_param), cudaMemcpyHostToDevice);
    //cudaMemcpy(data_gpu, data, (bits / 8)*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(sequenceGPU, sequence, MIN(MAXNUMOFTEMPLATES, numOfTemplates[m])*m*sizeof(double), cudaMemcpyHostToDevice);
    
    

    /////////////////////////////////////
    //线性复杂度计算所需要的参数
    /////////////////////////////////////
	N = 0;
	N = bits / M;
	int blocksPerGrid = (N + 512 -1) / 512;
	//申请host端内存及初始化
	double *h_partial_sum;	
    h_partial_sum = (double*)malloc( N*sizeof(double));

	time_t start;
	time_t end;
	double duration;
	gpu_param M_cpu;
	M_cpu.M = M;
	//int M[3] = {0};
	//M[0] = 500;
	int v[7] = {0};
	M_cpu.u = u;
	M_cpu.sign = sign;
	M_cpu.N = N;
//	printf("M:%d u:%f sign:%d\n",M_cpu.M,M_cpu.u,M_cpu.sign);
	
	
	//分配显存空间
	int size = sizeof(double);
	double *d_partial_sum;
	gpu_param *M_gpu;
	unsigned char *data_gpu;

//	start = clock();
	cudaMalloc((void**)&d_partial_sum,N*size);
	cudaMalloc((void**)&M_gpu,sizeof(gpu_param));
	cudaMalloc((void**)&data_gpu,(bits / 8)*sizeof(unsigned char));
	//printf("111111111\n");
	//把数据从Host传到Device
	cudaMemcpy(M_gpu, &M_cpu, sizeof(gpu_param), cudaMemcpyHostToDevice);
	cudaMemcpy(data_gpu, data, (bits / 8)*sizeof(unsigned char), cudaMemcpyHostToDevice);

	dim3 Grid(8,334,1);
	dim3 Block(148,6,1);
	//调用内核函数
	ReductionSum<<<Grid,Block>>>(WJ_gpu,CHI2_gpu,Mnon_gpu,sequenceGPU,d_partial_sum,M_gpu,data_gpu);
   // nonoverlapKernel<<<8,148>>>(WJ_gpu,CHI2_gpu,M_gpu,data_gpu,sequenceGPU);
    //nonoverlapKernel(int* WJ, double *CHI2_gpu, gpu_param *M_gpu, unsigned char* data, unsigned char*sequence)


	//将结果传回到主机端
	cudaMemcpy(WJ, WJ_gpu, 148*8*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_partial_sum, d_partial_sum, size*N, cudaMemcpyDeviceToHost);
    


	for(int i = 0;i<N;i++){
		if (h_partial_sum[i] <= -2.5)
		v[0]++;
	else if (h_partial_sum[i] <= -1.5)
		v[1]++;
	else if (h_partial_sum[i] <= -0.5)
		v[2]++;
	else if (h_partial_sum[i] <= 0.5)
		v[3]++;
	else if (h_partial_sum[i] <= 1.5)
		v[4]++;
	else if (h_partial_sum[i] <= 2.5)
		v[5]++;
	else
		v[6]++;

	}
	//end = clock();
	//printf("time:%f\n",(double)(end-start)/CLOCKS_PER_SEC);
		//Compute P-value
		double v_obs = 0.00;
	//	double p_value;
		double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
		for (int i = 0; i< 7; i++) {
			//printf("vi = %d\n",v[i]);

			v_obs += pow((double)v[i] - N*pi[i], 2) / (N*pi[i]);
		}
		//p_value = cephes_igamc(3.0, v_obs / 2.0);

  //	end = clock();
	  
	//duration = (double)(end-start)/CLOCKS_PER_SEC/ROUND;

	//printf("time1 = %f\n",duration);
	//printf("v_obs = %f\n",v_obs);
	int flag = 1;
	for( jj=0; jj<MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++ ) {
		sum = 0;
		chi2 = 0.0;                                   /* Compute Chi Square */
		//printf("wj:%d i:%d\n",WJ[jj],jj);
		for ( i=0; i<8; i++ ) {
		//	if ( m == 10 )
			//	fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
		//	else
			//	fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
			//printf("WJ:%d\n",WJ[i]);
			chi2 += pow(((double)WJ[jj*8+i] - lambda)/pow(varWj, 0.5), 2);
		}
		p_value = cephes_igamc(8/2.0, chi2/2.0);
		if (p_value < alpha) flag = 0;
	// 	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
	// 	printf("\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

	//	printf("%9.6f %f %s %3d\n", chi2, p_value, p_value < alpha ? "FAILURE" : "SUCCESS", jj);
	// 	if ( SKIP > 1 )
	// 		fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
	//}
		}
	
	//释放显存空间
	cudaFree(d_partial_sum);
	cudaFree(M_gpu);
	cudaFree(data_gpu);
 
    free(h_partial_sum);
    cudaDeviceReset();
	return v_obs;
    }

}