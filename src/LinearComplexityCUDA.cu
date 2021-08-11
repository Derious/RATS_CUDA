#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rats.h"

#define BITS	1000000
#define ROUND	1

typedef struct para{
	int M;
	double u;
	int sign;
	int N;
} gpu_param;

__device__ int GET_EPSILON(unsigned char*data,int offset){

	return ((data[offset/8]>>(7-(offset)%8))&1);
}

 
__global__ void ReductionSum(double *d_partial_sum, gpu_param *M_gpu, unsigned char* data)
{
	//申请共享内存，存在于每个block中 
	int M = M_gpu->M;
	double u = M_gpu->u;
	double Ti;
	int sign = M_gpu->sign;
	int N = M_gpu->N;
	//printf("M = %d\n",M);

	//__shared__ int v[7];
	//确定索引
	int threadId = blockIdx.x *blockDim.x + threadIdx.x;  
//	int i = threadIdx.x + blockIdx.x * blockDim.x;
	//int tid = threadIdx.x;
	int L , m , d , N_;

if(threadId<N){

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
			d = GET_EPSILON(data, threadId*M + N_);
			//printf("d: %d\n",d);
			for (int j = 1; j <= L; j++)
				d += C[j] * GET_EPSILON(data, threadId*M + N_ - j);
			if (d = d % 2)
			{
				for (int j = 0; j < M; j++)
				{
					T[j] = C[j];
					P[j] = 0;
				}
				for (int j = 0; j < M; j++)
				{
					if (B_[j] == 1)
						P[j + N_ - m] = 1;
				}
				for (int j = 0; j < M; j++)
					C[j] = (C[j] + P[j]) % 2;
				if (L <= N_ / 2)
				{
					L = N_ + 1 - L;
					m = N_;
					for (int j = 0; j<M; j++)
						B_[j] = T[j];
				}
			}
		}
		Ti = sign * (L - u) + 2.0 / 9.0;
		d_partial_sum[threadId] = Ti;
		// printf("Ti: %f\n",Ti);


	//printf("M:%d u:%f sign:%d Ti:%f\n",M,u,sign,Ti);
	free(B_);
	free(P);
	free(C);
	free(T);



}
}
 
double LinearComplexity(double alpha, unsigned char *data, int bits, int M, double u, int sign, LinearComplexity_V *value)
{
	int N = 0;
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


	//调用内核函数
	ReductionSum<<<blocksPerGrid,512>>>(d_partial_sum,M_gpu,data_gpu);


	//将结果传回到主机端
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
	
	//释放显存空间
	cudaFree(d_partial_sum);
	cudaFree(M_gpu);
	cudaFree(data_gpu);
 
    free(h_partial_sum);
	return v_obs;

}