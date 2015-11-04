#include "dc2.h"
#include "dc2-internal.h"
#include <cassert>
#include <cmath>
#include <cstring>
using namespace std;

template <>
void Sample2<FT_UNIT_SPHERE>(
    const Function& f,
    const float* x, const float* y, float* s, int n) {
    assert(f.mFunctionType == FT_UNIT_SPHERE);

    for (int i = 0; i < n; ++i) {
        s[i] = x[i] * x[i] + y[i] * y[i] - 1.0f;
    }
}

template <>
void SampleGradient2<FT_UNIT_SPHERE>(
    const Function& f,
    const float* x, const float* y,
    float* ds0, float* ds1,
    int n) {
    assert(f.mFunctionType == FT_UNIT_SPHERE);

    const float dh = 1e-5f;

    for (int i = 0; i < n; ++i) {
        ds0[i] = 2 * x[i];
        ds1[i] = 2 * y[i];
        float l = sqrt(ds0[i] * ds0[i] + ds1[i] * ds1[i]);
        ds0[i] /= l;
        ds1[i] /= l;
    }
}

template <>
void CollectIntersectionEdges2<FT_UNIT_SPHERE>(
    const Function& f,
    const float* x0, const float* y0,
    const float* x1, const float* y1,
    const float* v0, const float* v1, int n,
    float* xlow, float* ylow,
    float* xhigh, float* yhigh,
    int* ens,
    int* en) {
    assert(f.mFunctionType == FT_UNIT_SPHERE);

    int top = 0;
    for (int i = 0; i < n - 1; ++i) {
        int count = 0;
        int which = 0;

        if (v0[i] < 0.0f && v0[i + 1] >= 0.0f ||
            v0[i] >= 0.0f && v0[i + 1] < 0.0f) {
            const bool flag = v0[i] < 0.0f;
            xlow[top] = flag ? x0[i] : x0[i + 1];
            ylow[top] = flag ? y0[i] : y0[i + 1];
            xhigh[top] = flag ? x0[i + 1] : x0[i];
            yhigh[top] = flag ? y0[i + 1] : y0[i];
            which = 0;
            ++count;
            ++top;
        }

        if (v0[i] < 0.0f && v1[i] >= 0.0f ||
            v0[i] >= 0.0f && v1[i] < 0.0f) {
            const bool flag = v0[i] < 0.0f;
            xlow[top] = flag ? x0[i] : x1[i];
            ylow[top] = flag ? y0[i] : y1[i];
            xhigh[top] = flag ? x1[i] : x0[i];
            yhigh[top] = flag ? y1[i] : y0[i];
            which = 1;
            ++count;
            ++top;
        }

        ens[i] = 0;
        if (count == 2)
            ens[i] = 3; // 11b
        else if (count == 1)
            ens[i] = 1 << which;
    }
    *en = top;
}

template <>
void SolveIntersection2<FT_UNIT_SPHERE>(
    const Function& f,
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n) {
    assert(f.mFunctionType == FT_UNIT_SPHERE);

    const float eps = 1e-12f;

    for (int i = 0; i < n; ++i) {
        float xh = xhigh[i];
        float xl = xlow[i];
        float yh = yhigh[i];
        float yl = ylow[i];

        float cx = (xl + xh) * 0.5f;
        float cy = (yl + yh) * 0.5f;
        float v = cx * cx + cy * cy - 1.0f;

        while (abs(v) >= eps) {
            if (v > 0) {
                xh = cx;
                yh = cy;
            }
            else {
                xl = cx;
                yl = cx;
            }
            v = cx * cx + cy * cy - 1.0f;
        }
        x[i] = cx;
        y[i] = cy;
    }
}

void ConstructQEF2(
    const float* ix0, const float* iy0,
    const float* ix1, const float* iy1,
    const float* nx0, const float* ny0,
    const float* nx1, const float* ny1,
    const int* ens0, const int* ens1, int n,
    float* f, bool* h, int* m) {
    int count_total[] = {0, 1, 1, 2};
    int count_horizontal[] = {0, 1, 0, 1};
    int count_vertical[] = {0, 0, 1, 1};

    *m = 0;
    int start0 = 0, start1 = 0;
    for (int i = 0; i < n - 1; ++i) {
        int c0 = count_total[ens0[i]];
        int c1 = count_horizontal[ens1[i]];
        int c2 = count_total[ens1[i]];
        int c3 = count_vertical[ens0[i + 1]];

        if (c0 + c1 + c3 < 2){
            h[i] = false;
            start0 += c0;
            start1 += c2;
            continue;
        }

        h[i] = true;

        float A[8];
        float b[4];
        memset(A, 0, sizeof(A));
        memset(b, 0, sizeof(b));

        float gx = 0.0f;
        float gy = 0.0f;

        for (int j = 0; j < c0; ++j) {
            const float nx = nx0[start0 + j];
            const float ny = ny0[start0 + j];
            const float ix = ix0[start0 + j];
            const float iy = iy0[start0 + j];
            A[j * 2] = nx;
            A[j * 2 + 1] = ny;
            b[j] = nx * ix + ny * iy;
            gx += ix / (c0 + c1 + c3);
            gy += iy / (c0 + c1 + c3);
        }

        for (int j = 0; j < c1; ++j) {
            const float nx = nx1[start1 + j];
            const float ny = ny1[start1 + j];
            const float ix = ix1[start1 + j];
            const float iy = ix1[start1 + j];
            A[(c0 + j) * 2] = nx;
            A[(c0 + j) * 2 + 1] = ny;
            b[c0 + j] = nx * ix + ny * iy;
            gx += ix / (c0 + c1 + c3);
            gy += iy / (c0 + c1 + c3);
        }

        for (int j = 0; j < c3; ++j) {
            const float nx = nx0[start0 + c0 + j];
            const float ny = ny0[start0 + c0 + j];
            const float ix = ix0[start0 + c0 + j];
            const float iy = iy0[start0 + c0 + j];
            A[(c0 + c1 + j) * 2] = nx;
            A[(c0 + c1 + j) * 2 + 1] = ny;
            b[c0 + c1 + j] = mx * ix + ny * iy;
            gx += ix / (c0 + c1 + c3);
            gy += iy / (c0 + c1 + c3);
        }

        // Compute ATA and ATb.
        const float f0 = A[0] * A[0] + A[2] * A[2] + A[4] * A[4] + A[6] * A[6];
        const float f1 = A[0] * A[1] + A[2] * A[3] + A[4] * A[5] + A[6] * A[7];
        const float f2 = A[1] * A[1] + A[3] * A[3] + A[5] * A[5] + A[7] * A[7];

        const float f3 = A[0] * b[0] + A[2] * b[1] + A[4] * b[2] + A[6] * b[3];
        const float f4 = A[1] * b[0] + A[3] * b[1] + A[5] * b[2] + A[7] * b[3];

        float data[] = {
            f0, f1, f2, f3, f4, gx, gy
        };
        memcpy(f + i * 7, data, sizeof(data));
        *m += 1;

        start0 += c0;
        start1 += c2;
    }
}

void SolveQEF2(const float* f, float* p, int n) {

}
