
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_constants.h"
#include "math_functions.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

typedef struct Orb{
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double m;
	double st;
}Orb;

#define G 0.000005
#define MIN_DIST 1.0
const int unitSize = 8;//unitSize:每个计算单位的大小，多少个float
double wide = 10000;
double mass = 10;
double velo = 0.005;
int saveTimes = 0;

void calcGravity(Orb*o, Orb*ta, double dist, double*gx, double*gy, double*gz);
void calcOne(Orb*o, int oId, Orb*olist, int nUnit);
Orb* newOrbList(int nUnit, int style);
void deleteOrbList(Orb *olist);
void initOrbList(Orb *olist, int nUnit, int style);

/* calc gravity between two*/
void calcGravity(Orb*o, Orb*ta, double dist, double*gx, double*gy, double*gz) {
	double a = ta->m / (dist*dist) * G;
	*gx += -a * (o->x - ta->x) / dist;
	*gy += -a * (o->y - ta->y) / dist;
	*gz += -a * (o->z - ta->z) / dist;
}
/* Orb update once with list */
void calcOne(Orb*o, int oId, Orb*olist, int nUnit) {
	//if (o->st < 0) {
	int i = 0, isTooRappid = 0;
	double gax = 0, gay = 0, gaz = 0, dist = 0;
	for (i = 0; i<nUnit; ++i) {
		Orb* ta = olist + i;
		if (o->st < 0 && ta->st < 0 && oId != i) {
			dist = sqrt(((ta->x - o->x)*(ta->x - o->x) + (ta->y - o->y)*(ta->y - o->y) + (ta->z - o->z)*(ta->z - o->z)));
			isTooRappid = dist*dist<(o->vx*o->vx + o->vy*o->vy + o->vz*o->vz) * 10;
			if (dist<MIN_DIST || isTooRappid) {
				// crash
				if (o->m < ta->m) {
					o->st = -o->st;
					printf("one crash: oid=%d, tid=%d dist=%f isTooRappid=%d\n", oId, i, dist, isTooRappid);
				}
				continue;
			}
			calcGravity(o, ta, dist, &gax, &gay, &gaz);
		}
	}
	o->x += o->vx;
	o->y += o->vy;
	o->z += o->vz;
	o->vx += gax;
	o->vy += gay;
	o->vz += gaz;
	//}
}

Orb* newOrbList(int nUnit, int style) {

	int listSizeOfByte = sizeof(double) * unitSize * nUnit;
	Orb* list = (Orb*)malloc(listSizeOfByte);
	initOrbList(list, nUnit, style);
	return list;
}

void deleteOrbList(Orb *olist) {
	if (olist != NULL) {
		delete(olist);
	}
}
void initOrbList(Orb *olist, int nUnit, int style) {
	int i = 0;
	srand(time(NULL));
	for (i = 0; i<nUnit; ++i) {
		//Orb* o = (Orb*)(list+i*unitSize);
		Orb* o = olist + i;
		o->x = (double)rand() / (double)RAND_MAX*wide - wide / 2.0;
		o->y = (double)rand() / (double)RAND_MAX*wide - wide / 2.0;
		o->z = (double)rand() / (double)RAND_MAX*wide - wide / 2.0;
		o->vx = (double)rand() / (double)RAND_MAX*velo - velo / 2.0;
		o->vy = (double)rand() / (double)RAND_MAX*velo - velo / 2.0;
		o->vz = (double)rand() / (double)RAND_MAX*velo - velo / 2.0;
		o->m = (double)rand() / (double)RAND_MAX*mass;
		o->st = -(double)i;
	}
}

int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

cudaError_t calcOrbsWithCuda(Orb* olist, int nUnit, int nTimes) {
	int* dev_a = 0;
	int* dev_b = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_a, nUnit * sizeof(Orb));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

Error:
	cudaFree(dev_a);
	cudaFree(dev_b);

	return cudaStatus;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
