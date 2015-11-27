#include "cudc2.h"
#include <cstdlib>
#include <cstdio>
using namespace std;

#define TILE_WIDTH 32
#define TILE_HEIGHT 32
#define BLOCK_SIZE (TILE_WIDTH * TILE_HEIGHT)

template <FunctionType FT>
__device__ float Sample2(float x, float y);

template <>
__device__ float Sample2<FT_UNIT_SPHERE>(float x, float y) {
    return x * x + y * y - 1.0f;
}

__device__ bool IsIntersected(float a, float b) {
    return a < 0 && b >= 0 || a >= 0 && b < 0;
}

template <FunctionType FT>
__global__ void KDualContour2(
    float xs, float xt, int xn,
    float ys, float yt, int yn,
    float* p, int pcap, int* pcnt) {
    const float xstep = (xt - xs) / (xn - 1);
    const float ystep = (yt - ys) / (yn - 1);

    const int tx = blockIdx.x * TILE_WIDTH + threadIdx.x;
    const int ty = threadIdx.y;
    const int ti = threadIdx.y * TILE_WIDTH + threadIdx.x;
    const bool active = tx < (xn - 1) && ty < (yn - 1);

    const float x0 = xs + tx * xstep;
    const float x1 = x0 + xstep;
    const float y0 = ys + ty * ystep;
    const float y1 = y0 + ystep;

    const float s0 = Sample2<FT>(x0, y0);
    const float s1 = Sample2<FT>(x1, y0);
    const float s2 = Sample2<FT>(x0, y1);
    const float s3 = Sample2<FT>(x1, y1);

    int intersections = IsIntersected(s0, s1) +
        IsIntersected(s0, s2) +
        IsIntersected(s2, s3) +
        IsIntersected(s1, s3);

    int pred = active && intersections >= 2;
    const int total = __syncthreads_count(pred);

    __shared__ int base;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
        base = atomicAdd(pcnt, total);
    }

    __shared__ int s[BLOCK_SIZE];
    s[ti] = pred;

    __syncthreads();

    // Reduction step.
    for (int stride = 1; stride <= BLOCK_SIZE / 2; stride *= 2) {
        int index = (ti + 1) * stride * 2 - 1;
        if (index < BLOCK_SIZE)
            s[index] += s[index - stride];

        __syncthreads();
    }

    // Post scan step.
    for (int stride = BLOCK_SIZE / 4; stride > 0; stride /= 2) {
        int index = (ti + 1) * stride * 2 - 1;
        if (index + stride < BLOCK_SIZE)
            s[index + stride] += s[index];

        __syncthreads();
    }

    // Generate the vertex.
    if (pred) {
        const float cx = (x0 + x1) * 0.5f;
        const float cy = (y0 + y1) * 0.5f;
        p[(base + s[ti] - 1) * 2] = cx;
        p[(base + s[ti] - 1) * 2 + 1] = cy;
    }
}

bool DualContour2(
    FunctionType ft,
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap, int* pcnt,
    int* edges, int ecap, int* ecnt) {
    int blocks = (n - 1) / TILE_WIDTH;
    if ((n - 1) % TILE_WIDTH != 0)
        ++blocks;

    int tiles = (n - 1) / TILE_HEIGHT;
    if ((n - 1) % TILE_HEIGHT != 0)
        ++tiles;

    const float tilestep = (yt - ys) / (n - 1) * TILE_HEIGHT;

    dim3 blockDim(TILE_WIDTH, TILE_HEIGHT);
    dim3 gridDim(blocks, 1);

    // Allocate device memory.
    float* pd = nullptr;
    cudaMalloc(&pd, sizeof(float) * 2 * pcap);
    int* pcntd = nullptr;
    cudaMalloc(&pcntd, sizeof(int));
    cudaMemset(pcntd, 0, sizeof(int));

    // Launch computation grid tile by tile.
    if (ft == FT_UNIT_SPHERE) {
        for (int tileY = 0; tileY < tiles; ++tileY) {
            KDualContour2<FT_UNIT_SPHERE><<<gridDim, blockDim>>>(
                xs, xt, n,
                ys + tilestep * tileY, ys + tilestep * (tileY + 1),
                TILE_HEIGHT + 1,
                p, pcap, pcnt);

            cudaDeviceSynchronize();
        }
    }

    // Fetch device data.
    cudaMemcpy(pcnt, pcntd, sizeof(int),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(p, pd, sizeof(float) * 2 * pcap,
        cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    // Release resources.
    cudaFree(pcntd);
    cudaFree(pd);

    return true;
}