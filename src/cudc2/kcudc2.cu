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

    const float a = (Sample2<FT>(x + dh, y) - Sample2<FT>(x - dh, y)) * rev2dh;
    const float b = (Sample2<FT>(x, y + dh) - Sample2<FT>(x, y - dh)) * rev2dh;
    const float l = sqrt(a * a + b * b);

    *dx = a / l;
    *dy = b / l;
}

__device__ bool IsIntersected(float a, float b) {
    return a < 0 && b >= 0 || a >= 0 && b < 0;
}

template <FunctionType FT>
__device__ void SolveIntersection2(
    float xlow, float ylow,
    float xhigh, float yhigh,
    float* x, float* y) {
    const float EPS = 1e-6f;

    float mx = (xlow + xhigh) * 0.5f;
    float my = (ylow + yhigh) * 0.5f;
    float v = Sample2<FT>(mx, my);

    while (fabs(v) >= EPS) {
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
    *x = mx;
    *y = my;
}

__device__ void SVD2(
    const float* M,
    float* U,
    float* S,
    float* VT) {
    const float e = 0.5f * (M[0] + M[3]);
    const float f = 0.5f * (M[0] - M[3]);
    const float g = 0.5f * (M[1] + M[2]);
    const float h = 0.5f * (M[1] - M[2]);

    // one half of w1 + w2
    const float hw1pw2 = sqrt(e * e + h * h);
    const float hw1mw2 = sqrt(f * f + g * g);

    const float atangf = g == 0.0f && f == 0.0f ? 0.0f : atan2(g, f);
    const float atanhe = h == 0.0f && e == 0.0f ? 0.0f : atan2(h, e);
    
    S[0] = hw1pw2 + hw1mw2;
    S[1] = hw1pw2 - hw1mw2;

    const float b = (atanhe - atangf) * 0.5f;
    const float r = (atanhe + atangf) * 0.5f;

    const float cb = cos(b);
    const float sb = sin(b);
    const float cr = cos(r);
    const float sr = sin(r);

    U[0] = cb;
    U[1] = sb;
    U[2] = -sb;
    U[3] = cb;

    VT[0] = cr;
    VT[1] = sr;
    VT[2] = -sr;
    VT[3] = cr;
}

__device__ void PseudoInverse2(
    const float* M,
    float* pinv) {
    const float EPS = 1e-7f;

    float U[4], S[2], VT[4];
    SVD2(M, U, S, VT);

    S[0] = abs(S[0]) < EPS ? 0.0f : 1.0f / S[0];
    S[1] = abs(S[1]) < EPS ? 0.0f : 1.0f / S[1];

    const float v0 = S[0] * VT[0];
    const float v1 = S[1] * VT[2];
    const float v2 = S[0] * VT[1];
    const float v3 = S[1] * VT[3];

    pinv[0] = v0 * U[0] + v1 * U[1];
    pinv[1] = v0 * U[2] + v1 * U[3];
    pinv[2] = v2 * U[0] + v3 * U[1];
    pinv[3] = v2 * U[2] + v3 * U[3];
}

template <FunctionType FT>
__global__ void KDualContour2(
    float xs, float xt, int xn,
    float ys, float yt, int yn,
    float* p, int pcap, int* pcnt,
    int* ibuffer, int tileY) {
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

    // Compute the index of generated vertex.
    const int pi = base + s[ti] - 1;

    // Initialize the index buffer.
    const int which = tileY % 2;
    const int ibufferstride = (xn - 1) * TILE_HEIGHT;
    const int bi = which * ibufferstride + ty * (xn - 1) + tx;
    if (active)
        ibuffer[bi] = -1;

    // Generate the vertex.
    if (pred) {
        int k = 0;
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

            A[2 * k] = dx;
            A[2 * k + 1] = dy;

            cx += x;
            cy += y;

            b[k] = dx * x + dy * y;

            ++k;
        }

        if (IsIntersected(s0, s2)) {
            const bool n = s0 < s2;
            const float ylow = n ? y0 : y1;
            const float yhigh = n ? y1 : y0;
            float x, y;
            SolveIntersection2<FT>(x0, ylow, x0, yhigh, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * k] = dx;
            A[2 * k + 1] = dy;

            cx += x;
            cy += y;

            b[k] = dx * x + dy * y;

            ++k;
        }

        if (IsIntersected(s2, s3)) {
            const bool n = s2 < s3;
            const float xlow = n ? x0 : x1;
            const float xhigh = n ? x1 : x0;
            float x, y;
            SolveIntersection2<FT>(xlow, y1, xhigh, y1, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * k] = dx;
            A[2 * k + 1] = dy;

            cx += x;
            cy += y;

            b[k] = dx * x + dy * y;

            ++k;
        }

        if (IsIntersected(s1, s3)) {
            const bool n = s1 < s3;
            const float ylow = n ? y0 : y1;
            const float yhigh = n ? y1 : y0;
            float x, y;
            SolveIntersection2<FT>(x1, ylow, x1, yhigh, &x, &y);
            float dx, dy;
            SampleGradient2<FT>(x, y, &dx, &dy);

            A[2 * k] = dx;
            A[2 * k + 1] = dy;

            cx += x;
            cy += y;

            b[k] = dx * x + dy * y;

            ++k;
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
        f[5] = cx / k;
        f[6] = cy / k;

        // Solve QEF.
        float pinv[4];
        float ATA[4] = {f[0], f[1], f[1], f[2]};
        PseudoInverse2(ATA, pinv);
        const float ATAg[2] = {
            ATA[0] * f[5] + ATA[1] * f[6],
            ATA[2] * f[5] + ATA[3] * f[6]
        };
        const float d[2] = {
            f[3] - ATAg[0],
            f[4] - ATAg[1]
        };
        const float c[2] = {
            pinv[0] * d[0] + pinv[1] * d[1],
            pinv[2] * d[0] + pinv[3] * d[1]
        };
        const float px = c[0] + f[5];
        const float py = c[1] + f[6];

        p[2 * pi] = px;
        p[2 * pi + 1] = py;

        ibuffer[bi] = pi;
    }
}

__global__ void KBuildTopology(
    const int* ibuffer, int tileY,
    int xn, int yn,
    int* edges, int ecap, int* ecnt) {
    const int which = tileY % 2;
    const int ibufferstride = (xn - 1) * TILE_HEIGHT;
    const int* ib = ibuffer + which * ibufferstride;
    const int* prev = ibuffer + (1 - which) * ibufferstride;

    const int tx = blockIdx.x * TILE_WIDTH + threadIdx.x;
    const int ty = threadIdx.y;
    const int ti = threadIdx.y * TILE_WIDTH + threadIdx.x;
    const bool active = tx < (xn - 1) && ty < (yn - 1);

    const int self = active ? ib[ty * (xn - 1) + tx] : -1;
    const int left = tx == 0 ? -1 : ib[ty * (xn - 1) + tx - 1];
    const int down = (active && !(tileY == 0 && ty == 0)) ?
        (ty == 0 ? prev[(TILE_HEIGHT - 1) * (xn - 1) + tx] :
        ib[(ty - 1) * (xn - 1) + tx]) : -1;

    int count = 0;
    if (self != -1) {
        if (left != -1)
            ++count;
        if (down != -1)
            ++count;
    }

    const int total1 = __syncthreads_count(count == 1);
    const int total2 = __syncthreads_count(count == 2);
    const int total = total1 + total2 * 2;

    __shared__ int base;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
        base = atomicAdd(ecnt, total);
    }

    __shared__ int s[BLOCK_SIZE];
    s[ti] = count;

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

    const int ei = base + s[ti] - count;

    if (count != 0 && active) {
        const int node = left != -1 ? left : down;
        edges[2 * ei] = self;
        edges[2 * ei + 1] = node;
        --count;
    }
    if (count != 0 && active) {
        edges[2 * (ei + 1)] = self;
        edges[2 * (ei + 1) + 1] = down;
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
    int r = (n - 1) % TILE_HEIGHT;
    if (r != 0)
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
    
    int* edgesd = nullptr;
    cudaMalloc(&edgesd, sizeof(int) * 2 * ecap);
    int* ecntd = nullptr;
    cudaMalloc(&ecntd, sizeof(int));
    cudaMemset(ecntd, 0, sizeof(int));

    // Allocate double index buffer.
    const int ibuffersize = sizeof(int) * (n - 1) * TILE_HEIGHT * 2;
    int* ibufferd = nullptr;
    cudaMalloc(&ibufferd, ibuffersize);
    cudaMemset(ibufferd, -1, ibuffersize);

    // Launch computation grid tile by tile.
    if (ft == FT_UNIT_SPHERE) {
        for (int tileY = 0; tileY < tiles; ++tileY) {
            int yn = TILE_HEIGHT + 1;
            if (tileY == tiles - 1)
                yn = r + 1;

            KDualContour2<FT_UNIT_SPHERE><<<gridDim, blockDim>>>(
                xs, xt, n,
                ys + tilestep * tileY, ys + tilestep * (tileY + 1),
                yn,
                pd, pcap, pcntd,
                ibufferd, tileY);

            cudaDeviceSynchronize();

            KBuildTopology<<<gridDim, blockDim>>>(
                ibufferd,
                tileY,
                n, yn,
                edgesd, ecap, ecntd);

            cudaDeviceSynchronize();
        }
    }

    // Fetch device data.
    cudaMemcpy(pcnt, pcntd, sizeof(int),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(p, pd, sizeof(float) * 2 * pcap,
        cudaMemcpyDeviceToHost);
    
    cudaMemcpy(ecnt, ecntd, sizeof(int),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(edges, edgesd, sizeof(int) * 2 * ecap,
        cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    // Release resources.
    cudaFree(ibufferd);
    cudaFree(ecntd);
    cudaFree(edgesd);
    cudaFree(pcntd);
    cudaFree(pd);

    return true;
}