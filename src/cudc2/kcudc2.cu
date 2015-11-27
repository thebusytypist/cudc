#include "cudc2.h"
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

template <FunctionType FT>
__device__ void SampleGradient2(float x, float y, float* dx, float* dy) {
    const float dh = 1e-5f;
    const float rev2dh = 0.5f / dh;

    *dx = (Sample2<FT>(x + dh, y) - Sample2<FT>(x - dh, y)) * rev2dh;
    *dy = (Sample2<FT>(x, y + dh) - Sample2<FT>(x, y - dh)) * rev2dh;
}

__device__ bool IsIntersected(float a, float b) {
    return a < 0 && b >= 0 || a >= 0 && b < 0;
}

template <FunctionType FT>
__device__ void SolveIntersection2(
    float xlow, float ylow,
    float xhigh, float yhigh,
    float* x, float* y) {
    const float eps = 1e-6f;

    float mx = (xlow + xhigh) * 0.5f;
    float my = (ylow + yhigh) * 0.5f;
    float v = Sample2<FT>(mx, my);

    while (fabs(v) >= eps) {
        if (v > 0) {
            xhigh = mx;
            yhigh = my;
        }
        else {
            xlow = mx;
            ylow = my;
        }
        mx = (xlow + xhigh) * 0.5f;
        my = (ylow + yhigh) * 0.5f;
        v = Sample2<FT>(mx, my);
    }
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
        int d = 0;
        float cx = 0.0f, cy = 0.0f;
        float A[8]= {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        float b[4];

        if (IsIntersected(s0, s1)) {
            const bool n = s0 < s1;
            const float xlow = n ? x0 : x1;
            const float xhigh = n ? x1 : x0;
            float x, y;
            SolveIntersection2<FT>(xlow, y0, xhigh, y0, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * d] = dx;
            A[2 * d + 1] = dy;

            cx += x;
            cy += y;

            b[d] = dx * x + dy * y;

            ++d;
        }

        if (IsIntersected(s0, s2)) {
            const bool n = s0 < s2;
            const float ylow = n ? y0 : y1;
            const float yhigh = n ? y1 : y0;
            float x, y;
            SolveIntersection2<FT>(x0, ylow, x0, yhigh, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * d] = dx;
            A[2 * d + 1] = dy;

            cx += x;
            cy += y;

            b[d] = dx * x + dy * y;

            ++d;
        }

        if (IsIntersected(s2, s3)) {
            const bool n = s2 < s3;
            const float xlow = n ? x0 : x1;
            const float xhigh = n ? x1 : x0;
            float x, y;
            SolveIntersection2<FT>(xlow, y1, xhigh, y1, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * d] = dx;
            A[2 * d + 1] = dy;

            cx += x;
            cy += y;

            b[d] = dx * x + dy * y;

            ++d;
        }

        if (IsIntersected(s1, s3)) {
            const bool n = s1 < s3;
            const float ylow = n ? y0 : y1;
            const float yhigh = n ? y1 : y0;
            float x, y;
            SolveIntersection2<FT>(x1, ylow, x1, yhigh, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * d] = dx;
            A[2 * d + 1] = dy;

            cx += x;
            cy += y;

            b[d] = dx * x + dy * y;

            ++d;
        }

        // 7 floats for QEF equation.
        float f[7];

        // Compute ATA.
        f[0] = A[0] * A[0] + A[2] * A[2] + A[4] * A[4] + A[6] * A[6];
        f[1] = A[0] * A[1] + A[2] * A[3] + A[4] * A[5] + A[6] * A[7];
        f[2] = A[1] * A[1] + A[3] * A[3] + A[5] * A[5] + A[7] * A[7];

        // Compute ATb.
        f[3] = A[0] * b[0] + A[2] * b[1] + A[4] * b[2] + A[6] * b[3];
        f[4] = A[1] * b[0] + A[3] * b[1] + A[5] * b[2] + A[7] * b[3];

        // Compute mass point.
        f[5] = cx / d;
        f[6] = cy / d;

        // Solve QEF.
        
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
                pd, pcap, pcntd);

            cudaDeviceSynchronize();
        }
    }

    // Fetch device data.
    cudaMemcpy(pcnt, pcntd, sizeof(int),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(p, pd, sizeof(float) * 2 * pcap,
        cudaMemcpyDeviceToHost);
    *ecnt = 0;

    cudaDeviceSynchronize();

    // Release resources.
    cudaFree(pcntd);
    cudaFree(pd);

    return true;
}