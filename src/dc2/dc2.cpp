#include "dc2.h"
#include "dc2-internal.h"
#include <cmath>
#include <cassert>
using namespace std;

template <>
void Sample2<FT_UNIT_SPHERE>(
    const Function& f,
    const float* x, const float* y, float* s, int n) {
    assert(f.FunctionType == FT_UNIT_SPHERE);

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
    assert(f.FunctionType == FT_UNIT_SPHERE);

    const float dh = 1e-5f;

    for (int i = 0; i < n; ++i) {
        ds0[i] = 2 * x[i];
        ds1[i] = 2 * y[i];
    }
}

template <>
void CollectIntersectionEdges<FT_UNIT_SPHERE>(
    const Function& f,
    const float* x0, const float* y0,
    const float* x1, const float* y1,
    const float* v0, const float* v1, int n,
    float* xlow, float* ylow,
    float* xhigh, float* yhigh,
    int* ens,
    int* en) {
    assert(f.FunctionType == FT_UNIT_SPHERE);

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
void SolveIntersection<FT_UNIT_SPHERE>(
    const Function& f,
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n) {
    assert(f.FunctionType == FT_UNIT_SPHERE);

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
