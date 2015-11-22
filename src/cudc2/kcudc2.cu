#include <cstdlib>
#include <cstdio>
using namespace std;

#define BLOCK_SIZE 32

__global__ void StreamingCompactionKernel(
    int* blockOffset,
    float* data, int n,
    float* output, int cap) {
    const int i = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    int pred = data[i] < 0.5f;
    const int count = __syncthreads_count(pred);

    __shared__ int base;

    if (threadIdx.x == 0) {
        base = atomicAdd(blockOffset, count);
    }

    __shared__ int s[BLOCK_SIZE];
    s[threadIdx.x] = pred;

    __syncthreads();

    for (int stride = 1; stride <= BLOCK_SIZE / 2; stride *= 2) {
        int index = (threadIdx.x + 1) * stride * 2 - 1;
        if (index < BLOCK_SIZE)
            s[index] += s[index - stride];

        __syncthreads();
    }

    for (int stride = BLOCK_SIZE / 4; stride > 0; stride /= 2) {
        int index = (threadIdx.x + 1) * stride * 2 - 1;
        if (index + stride < BLOCK_SIZE)
            s[index + stride] += s[index];

        __syncthreads();
    }

    if (pred) {
        output[base + s[threadIdx.x] - 1] = data[i];
    }
}

bool StreamingCompact() {
    const int n = 64;
    const int capacity = n;

    float* data = (float*)malloc(sizeof(float) * n);
    for (int i = 0; i < n; ++i) {
        data[i] = (float)rand() / (float)RAND_MAX;
        printf("%f ", data[i]);
    }
    printf("\n\n");

    float* dataDevice = nullptr;
    cudaMalloc(&dataDevice, sizeof(float) * n);
    cudaMemcpy(dataDevice, data, sizeof(float) * n, cudaMemcpyHostToDevice);

    int* counterDevice = nullptr;
    cudaMalloc(&counterDevice, sizeof(int));
    cudaMemset(counterDevice, 0, sizeof(int));

    float* outputDevice = nullptr;
    cudaMalloc(&outputDevice, sizeof(float) * capacity);
    cudaMemset(outputDevice, 0, sizeof(float) * capacity);

    int grids = n / BLOCK_SIZE;
    if (grids % BLOCK_SIZE != 0)
        ++grids;
    StreamingCompactionKernel<<<grids, BLOCK_SIZE>>>(
        counterDevice,
        dataDevice, n,
        outputDevice, capacity);

    cudaDeviceSynchronize();

    float* result = (float*)malloc(sizeof(float) * capacity);
    cudaMemcpy(result, outputDevice, sizeof(float) * capacity,
        cudaMemcpyDeviceToHost);

    for (int i = 0; i < capacity; ++i) {
        printf("%f ", result[i]);
    }
    printf("\n");

    cudaFree(outputDevice);
    cudaFree(counterDevice);
    cudaFree(dataDevice);
    free(result);
    free(data);

    return true;
}
