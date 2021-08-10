#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define BITS	1000000
#define ROUND	1
#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

typedef struct para{
	int M;
	double u;
	int sign;
	int N;
} gpu_param;

const int threadsPerBlock=512 ; 
const int N=2000;
const int blocksPerGrid = (N + threadsPerBlock -1) / threadsPerBlock;
 
__global__ void ReductionSum(float *d_a, float *d_partial_sum, gpu_param *M_gpu, unsigned char* data)
{
	//申请共享内存，存在于每个block中 
	int M = M_gpu->M;
	double u = M_gpu->u;
	double Ti;
	int sign = M_gpu->sign;
	int N = M_gpu->N;
	//printf("M = %d\n",M);
	__shared__ float partialSum[threadsPerBlock];
	//__shared__ int v[7];
	//确定索引
	int threadId = blockIdx.x *blockDim.x + threadIdx.x;  
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int tid = threadIdx.x;
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


	//printf("M:%d u:%f sign:%d Ti:%f\n",M,u,sign,Ti);
	free(B_);
	free(P);
	free(C);
	free(T);



}

	// }
	// //传global memory数据到shared memory
	// partialSum[tid]=d_a[i];



	// //传输同步
	// __syncthreads();
	
	// //在共享存储器中进行规约
	// for(int stride = blockDim.x/2; stride > 0; stride/=2)
	// {
	// 	if(tid<stride) partialSum[tid]+=partialSum[tid+stride];
	// 	__syncthreads();
	// }
	// //将当前block的计算结果写回输出数组
	// if(tid==0)  
	// 	d_partial_sum[blockIdx.x] = partialSum[0];
}
 
int main()
{
	////////////////////////////读写文件测试
	unsigned char data[BITS / 8];
	unsigned char temp[1024];
	int b, i, rlen;
	FILE *fp;
	int flag;
	clock_t time_s, time_e;
	printf("1111111111111\n");
	if (!(fp = fopen("./data.e", "rb")))
	{
		printf("File open error!\n");
	}
	printf("1111111111111\n");
	memset(data, 0, BITS / 8);
	b = 0;
	while (b < BITS)
	{
		rlen = fread(temp, 1, 1024, fp);
		for (i = 0; i < rlen && b < BITS; i++)
		{
			if (temp[i] == '0' || temp[i] == '1')
			{
				data[b / 8] = (data[b / 8] << 1) + (temp[i] - '0');
				b++;
			}
		}
	}

	//for(int i=0;i<10;i++) printf("data:%x\n",data[i]);


	///////////////////////////





	//申请host端内存及初始化
	double   *h_a,*h_partial_sum;	
    h_a = (double*)malloc( N*sizeof(double) );
    h_partial_sum = (double*)malloc( N*sizeof(double));

	int sign,u;
	time_t start;
	time_t end;
	double duration;
	gpu_param M_cpu;
	M_cpu.M = 500;
	int M[3] = {0};
	M[0] = 500;
	int v[7] = {0};
	M_cpu.sign = ((M_cpu.M + 1) & 1) ? -1 : 1;
	M_cpu.u = M_cpu.M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M_cpu.M) * (M_cpu.M / 3.0 + 2.0 / 9.0);
	M_cpu.sign = (M_cpu.M & 1) ? -1 : 1;
	M_cpu.N = 2000;
	//printf("M:%d u:%f sign:%d\n",M_cpu.M,M_cpu.u,M_cpu.sign);
	
	for (int i=0; i < N; ++i)  h_a[i] = 1;
	
	//分配显存空间
	int size = sizeof(double);
	float *d_a;
	float *d_partial_sum;
	gpu_param *M_gpu;
	unsigned char *data_gpu;
	start = clock();

	cudaError_t err;






	cudaMalloc((void**)&d_a,N*size);
	cudaMalloc((void**)&d_partial_sum,N*size);
	cudaMalloc((void**)&M_gpu,sizeof(gpu_param));
	cudaMalloc((void**)&data_gpu,(BITS / 8)*sizeof(unsigned char));
	printf("111111111\n");
	//把数据从Host传到Device
	cudaMemcpy(d_a, h_a, size*N, cudaMemcpyHostToDevice);
	cudaMemcpy(M_gpu, &M_cpu, sizeof(gpu_param), cudaMemcpyHostToDevice);
	cudaMemcpy(data_gpu, data, (BITS / 8)*sizeof(unsigned char), cudaMemcpyHostToDevice);


	//cudaMemcpyToSymbol(v, v_cpu, 7*sizeof(int));
	printf("111111111\n");

for(int i=0;i<ROUND;i++){
	//调用内核函数
	ReductionSum<<<5,512>>>(d_a,d_partial_sum,M_gpu,data_gpu);
}

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
		double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
		for (int i = 0; i< 7; i++)
			v_obs += pow((double)v[i] - N*pi[i], 2) / (N*pi[i]);
	//	p_value = cephes_igamc(3.0, v_obs / 2.0);

  	end = clock();
	  
	duration = (double)(end-start)/CLOCKS_PER_SEC/ROUND;

	printf("time1 = %f\n",duration);


	// //将部分和求和
	// long sum=0;
    // for (int i=0; i < blocksPerGrid; ++i)  sum += h_partial_sum[i];

	// cout<<"sum1="<<sum<<" time1="<<duration<<endl;

	// start = clock();
	// sum = 0;
	// for (long i = 0; i < N; i++)
	// {
	// 	sum+=h_a[i];
	// 	/* code */
	// }
	// end = clock();
	// duration = (double)(end-start)/CLOCKS_PER_SEC;
 	// cout<<"sum2="<<sum<<" time2="<<duration<<endl;
	
	//释放显存空间
	cudaFree(d_a);
	cudaFree(d_partial_sum);
 
	free(h_a);
    free(h_partial_sum);
 END:
    return 0;
}