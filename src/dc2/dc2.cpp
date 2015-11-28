#include "dc2.h"
#include "dc2-internal.h"
#include "pinv.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace std;

template <>
void Sample2<FT_UNIT_SPHERE>(
    float xs, float xt,
    float ys, float yt,
    float* s, int n) {
    const float xstep = (xt - xs) / (n - 1);
    const float ystep = (yt - ys) / (n - 1);
    for (int i = 0; i < n; ++i) {
        const float x = xs + i * xstep;
        const float y = ys + i * ystep;
        s[i] = x * x + y * y - 1.0f;
    }
}

template <>
void Sample2<FT_UNIT_SPHERE>(
    const float* x,
    const float* y,
    float dx,
    float dy,
    float* s,
    int n) {
    for (int i = 0; i < n; ++i) {
        const float xi = x[i] + dx;
        const float yi = y[i] + dy;
        s[i] = xi * xi + yi * yi - 1.0f;
    }
}

template <>
void Sample2<FT_HEART>(
    float xs, float xt,
    float ys, float yt,
    float* s, int n) {
    const float xstep = (xt - xs) / (n - 1);
    const float ystep = (yt - ys) / (n - 1);
    for (int i = 0; i < n; ++i) {
        const float x = xs + i * xstep;
        const float y = ys + i * ystep;
        const float x2 = x * x;
        const float y2 = y * y;
        const float t = x2 + y2 - 1.0f;
        s[i] = t * t * t - x2 * y2 * y;
    }
}

template <>
void Sample2<FT_HEART>(
    const float* x,
    const float* y,
    float dx,
    float dy,
    float* s,
    int n) {
    for (int i = 0; i < n; ++i) {
        const float xi = x[i] + dx;
        const float yi = y[i] + dy;
        const float x2 = xi * xi;
        const float y2 = yi * yi;
        const float t = x2 + y2 - 1.0f;
        s[i] = t * t * t - x2 * y2 * yi;
    }
}

template <FunctionType FT>
void SampleGradient2(
    const float* x, const float* y,
    float* ds0, float* ds1,
    int n) {
    const float dh = 1e-5f;
    const float rev2dh = 0.5f / dh;

    int n4 = n / 4;
    int m = n4 << 2;
    int rn = n - m;

    float xreg[4], yreg[4];
    float sreg0[4], sreg1[4];
    float lreg[4];
    for (int i = 0; i < m; i += 4) {
        xreg[0] = x[i];
        xreg[1] = x[i + 1];
        xreg[2] = x[i + 2];
        xreg[3] = x[i + 3];

        yreg[0] = y[i];
        yreg[1] = y[i + 1];
        yreg[2] = y[i + 2];
        yreg[3] = y[i + 3];

        // Compute df/dx.
        Sample2<FT>(xreg, yreg, dh, 0.0f, sreg1, 4);
        Sample2<FT>(xreg, yreg, -dh, 0.0f, sreg0, 4);
        ds0[i] = (sreg1[0] - sreg0[0]) * rev2dh;
        ds0[i + 1] = (sreg1[1] - sreg0[1]) * rev2dh;
        ds0[i + 2] = (sreg1[2] - sreg0[2]) * rev2dh;
        ds0[i + 3] = (sreg1[3] - sreg0[3]) * rev2dh;

        // Compute df/dy.
        Sample2<FT>(xreg, yreg, 0.0f, dh, sreg1, 4);
        Sample2<FT>(xreg, yreg, 0.0f, -dh, sreg0, 4);
        ds1[i] = (sreg1[0] - sreg0[0]) * rev2dh;
        ds1[i + 1] = (sreg1[1] - sreg0[1]) * rev2dh;
        ds1[i + 2] = (sreg1[2] - sreg0[2]) * rev2dh;
        ds1[i + 3] = (sreg1[3] - sreg0[3]) * rev2dh;

        // Normalize the gradient vector.
        lreg[0] = sqrt(ds0[i] * ds0[i] + ds1[i] * ds1[i]);
        lreg[1] = sqrt(ds0[i + 1] * ds0[i + 1] + ds1[i + 1] * ds1[i + 1]);
        lreg[2] = sqrt(ds0[i + 2] * ds0[i + 2] + ds1[i + 2] * ds1[i + 2]);
        lreg[3] = sqrt(ds0[i + 3] * ds0[i + 3] + ds1[i + 3] * ds1[i + 3]);

        ds0[i] /= lreg[0]; ds1[i] /= lreg[0];
        ds0[i + 1] /= lreg[1]; ds1[i + 1] /= lreg[1];
        ds0[i + 2] /= lreg[2]; ds1[i + 2] /= lreg[2];
        ds0[i + 3] /= lreg[3]; ds1[i + 3] /= lreg[3];
    }

    // Handle the boundary case.
    xreg[0] = 0 < rn ? x[m] : 0.0f;
    xreg[1] = 1 < rn ? x[m + 1] : 0.0f;
    xreg[2] = 2 < rn ? x[m + 2] : 0.0f;
    xreg[3] = 0.0f;

    yreg[0] = 0 < rn ? y[m] : 0.0f;
    yreg[1] = 1 < rn ? y[m + 1] : 0.0f;
    yreg[2] = 2 < rn ? y[m + 2] : 0.0f;
    yreg[3] = 0.0f;

    // Compute df/dx.
    Sample2<FT>(xreg, yreg, dh, 0.0f, sreg1, 4);
    Sample2<FT>(xreg, yreg, -dh, 0.0f, sreg0, 4);
    if (0 < rn)
        ds0[m] = (sreg1[0] - sreg0[0]) * rev2dh;
    if (1 < rn)
        ds0[m + 1] = (sreg1[1] - sreg0[1]) * rev2dh;
    if (2 < rn)
        ds0[m + 2] = (sreg1[2] - sreg0[2]) * rev2dh;

    // Compute df/dy.
    Sample2<FT>(xreg, yreg, 0.0f, dh, sreg1, 4);
    Sample2<FT>(xreg, yreg, 0.0f, -dh, sreg0, 4);
    if (0 < rn)
        ds1[m] = (sreg1[0] - sreg0[0]) * rev2dh;
    if (1 < rn)
        ds1[m + 1] = (sreg1[1] - sreg0[1]) * rev2dh;
    if (2 < rn)
        ds1[m + 2] = (sreg1[2] - sreg0[2]) * rev2dh;

    // Normalize the gradient vector.
    lreg[0] = sqrt(ds0[m] * ds0[m] + ds1[m] * ds1[m]);
    lreg[1] = sqrt(ds0[m + 1] * ds0[m + 1] + ds1[m + 1] * ds1[m + 1]);
    lreg[2] = sqrt(ds0[m + 2] * ds0[m + 2] + ds1[m + 2] * ds1[m + 2]);
    lreg[3] = sqrt(ds0[m + 3] * ds0[m + 3] + ds1[m + 3] * ds1[m + 3]);

    if (0 < rn)
        ds0[m] /= lreg[0]; ds1[m] /= lreg[0];
    if (1 < rn)
        ds0[m + 1] /= lreg[1]; ds1[m + 1] /= lreg[1];
    if (2 < rn)
        ds0[m + 2] /= lreg[2]; ds1[m + 2] /= lreg[2];
}

void CollectIntersectionEdges2(
    float x0s, float x0t,
    float y0s, float y0t,
    float x1s, float x1t,
    float y1s, float y1t,
    const float* v0, const float* v1, int n,
    float* xlow, float* ylow,
    float* xhigh, float* yhigh,
    int* ens,
    int* en) {

    const float xstep = (x0t - x0s) / (n - 1);
    const float ystep = (y0t - y0s) / (n - 1);

    int top = 0;
    for (int i = 0; i < n - 1; ++i) {
        int count = 0;
        int which = 0;

        const float x0i = x0s + i * xstep;
        const float x0i1 = x0s + (i + 1) * xstep;
        const float y0i = y0s + i * ystep;
        const float y0i1 = y0s + (i + 1) * ystep;
        const float x1i = x1s + i * xstep;
        const float y1i = y1s + i * ystep;

        if (v0[i] < 0.0f && v0[i + 1] >= 0.0f ||
            v0[i] >= 0.0f && v0[i + 1] < 0.0f) {
            const bool flag = v0[i] < 0.0f;
            xlow[top] = flag ? x0i : x0i1;
            ylow[top] = flag ? y0i : y0i1;
            xhigh[top] = flag ? x0i1 : x0i;
            yhigh[top] = flag ? y0i1 : y0i;
            which = 0;
            ++count;
            ++top;
        }

        if (v0[i] < 0.0f && v1[i] >= 0.0f ||
            v0[i] >= 0.0f && v1[i] < 0.0f) {
            const bool flag = v0[i] < 0.0f;
            xlow[top] = flag ? x0i : x1i;
            ylow[top] = flag ? y0i : y1i;
            xhigh[top] = flag ? x1i : x0i;
            yhigh[top] = flag ? y1i : y0i;
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

template <FunctionType FT>
void SolveIntersection2(
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n) {
    const float eps = 1e-6f;

    int n4 = n / 4;
    int m = n4 << 2;
    int rn = n - m;

    float xlreg[4], xhreg[4];
    float ylreg[4], yhreg[4];
    float cxreg[4], cyreg[4];
    float vreg[4];
    int cmp[4];
    for (int i = 0; i < m; i += 4) {
        xlreg[0] = xlow[i];
        xlreg[1] = xlow[i + 1];
        xlreg[2] = xlow[i + 2];
        xlreg[3] = xlow[i + 3];

        xhreg[0] = xhigh[i];
        xhreg[1] = xhigh[i + 1];
        xhreg[2] = xhigh[i + 2];
        xhreg[3] = xhigh[i + 3];

        ylreg[0] = ylow[i];
        ylreg[1] = ylow[i + 1];
        ylreg[2] = ylow[i + 2];
        ylreg[3] = ylow[i + 3];

        yhreg[0] = yhigh[i];
        yhreg[1] = yhigh[i + 1];
        yhreg[2] = yhigh[i + 2];
        yhreg[3] = yhigh[i + 3];

        cxreg[0] = (xlreg[0] + xhreg[0]) * 0.5f;
        cxreg[1] = (xlreg[1] + xhreg[1]) * 0.5f;
        cxreg[2] = (xlreg[2] + xhreg[2]) * 0.5f;
        cxreg[3] = (xlreg[3] + xhreg[3]) * 0.5f;

        cyreg[0] = (ylreg[0] + yhreg[0]) * 0.5f;
        cyreg[1] = (ylreg[1] + yhreg[1]) * 0.5f;
        cyreg[2] = (ylreg[2] + yhreg[2]) * 0.5f;
        cyreg[3] = (ylreg[3] + yhreg[3]) * 0.5f;

        Sample2<FT>(cxreg, cyreg, 0.0f, 0.0f, vreg, 4);

        while (
            abs(vreg[0]) >= eps ||
            abs(vreg[1]) >= eps ||
            abs(vreg[2]) >= eps ||
            abs(vreg[3]) >= eps) {
            cmp[0] = vreg[0] > 0.0f;
            cmp[1] = vreg[1] > 0.0f;
            cmp[2] = vreg[2] > 0.0f;
            cmp[3] = vreg[3] > 0.0f;

            if (cmp[0]) {
                xhreg[0] = cxreg[0];
                yhreg[0] = cyreg[0];
            }
            else {
                xlreg[0] = cxreg[0];
                ylreg[0] = cyreg[0];
            }

            if (cmp[1]) {
                xhreg[1] = cxreg[1];
                yhreg[1] = cyreg[1];
            }
            else {
                xlreg[1] = cxreg[1];
                ylreg[1] = cyreg[1];
            }

            if (cmp[2]) {
                xhreg[2] = cxreg[2];
                yhreg[2] = cyreg[2];
            }
            else {
                xlreg[2] = cxreg[2];
                ylreg[2] = cyreg[2];
            }

            if (cmp[3]) {
                xhreg[3] = cxreg[3];
                yhreg[3] = cyreg[3];
            }
            else {
                xlreg[3] = cxreg[3];
                ylreg[3] = cyreg[3];
            }

            cxreg[0] = (xlreg[0] + xhreg[0]) * 0.5f;
            cxreg[1] = (xlreg[1] + xhreg[1]) * 0.5f;
            cxreg[2] = (xlreg[2] + xhreg[2]) * 0.5f;
            cxreg[3] = (xlreg[3] + xhreg[3]) * 0.5f;

            cyreg[0] = (ylreg[0] + yhreg[0]) * 0.5f;
            cyreg[1] = (ylreg[1] + yhreg[1]) * 0.5f;
            cyreg[2] = (ylreg[2] + yhreg[2]) * 0.5f;
            cyreg[3] = (ylreg[3] + yhreg[3]) * 0.5f;

            Sample2<FT>(cxreg, cyreg, 0.0f, 0.0f, vreg, 4);
        } // End of binary search.

        x[i] = cxreg[0];
        x[i + 1] = cxreg[1];
        x[i + 2] = cxreg[2];
        x[i + 3] = cxreg[3];

        y[i] = cyreg[0];
        y[i + 1] = cyreg[1];
        y[i + 2] = cyreg[2];
        y[i + 3] = cyreg[3];
    } // End of 4x cases.

    // Handle the boundary cases.
    xlreg[0] = 0 < rn ? xlow[m] : 0.0f;
    xlreg[1] = 1 < rn ? xlow[m + 1] : 0.0f;
    xlreg[2] = 2 < rn ? xlow[m + 2] : 0.0f;
    xlreg[3] = 0.0f;

    xhreg[0] = 0 < rn ? xhigh[m] : 0.0f;
    xhreg[1] = 1 < rn ? xhigh[m + 1] : 0.0f;
    xhreg[2] = 2 < rn ? xhigh[m + 2] : 0.0f;
    xhreg[3] = 0.0f;

    ylreg[0] = 0 < rn ? ylow[m] : 0.0f;
    ylreg[1] = 1 < rn ? ylow[m + 1] : 0.0f;
    ylreg[2] = 2 < rn ? ylow[m + 2] : 0.0f;
    ylreg[3] = 0.0f;

    yhreg[0] = 0 < rn ? yhigh[m] : 0.0f;
    yhreg[1] = 1 < rn ? yhigh[m + 1] : 0.0f;
    yhreg[2] = 2 < rn ? yhigh[m + 2] : 0.0f;
    yhreg[3] = 0.0f;

    cxreg[0] = (xlreg[0] + xhreg[0]) * 0.5f;
    cxreg[1] = (xlreg[1] + xhreg[1]) * 0.5f;
    cxreg[2] = (xlreg[2] + xhreg[2]) * 0.5f;
    cxreg[3] = (xlreg[3] + xhreg[3]) * 0.5f;

    cyreg[0] = (ylreg[0] + yhreg[0]) * 0.5f;
    cyreg[1] = (ylreg[1] + yhreg[1]) * 0.5f;
    cyreg[2] = (ylreg[2] + yhreg[2]) * 0.5f;
    cyreg[3] = (ylreg[3] + yhreg[3]) * 0.5f;

    Sample2<FT>(cxreg, cyreg, 0.0f, 0.0f, vreg, 4);

    while (
        0 < rn && abs(vreg[0]) >= eps ||
        1 < rn && abs(vreg[1]) >= eps ||
        2 < rn && abs(vreg[2]) >= eps) {
        cmp[0] = vreg[0] > 0.0f;
        cmp[1] = vreg[1] > 0.0f;
        cmp[2] = vreg[2] > 0.0f;
        cmp[3] = vreg[3] > 0.0f;

        if (cmp[0]) {
            xhreg[0] = cxreg[0];
            yhreg[0] = cyreg[0];
        }
        else {
            xlreg[0] = cxreg[0];
            ylreg[0] = cyreg[0];
        }

        if (cmp[1]) {
            xhreg[1] = cxreg[1];
            yhreg[1] = cyreg[1];
        }
        else {
            xlreg[1] = cxreg[1];
            ylreg[1] = cyreg[1];
        }

        if (cmp[2]) {
            xhreg[2] = cxreg[2];
            yhreg[2] = cyreg[2];
        }
        else {
            xlreg[2] = cxreg[2];
            ylreg[2] = cyreg[2];
        }

        cxreg[0] = (xlreg[0] + xhreg[0]) * 0.5f;
        cxreg[1] = (xlreg[1] + xhreg[1]) * 0.5f;
        cxreg[2] = (xlreg[2] + xhreg[2]) * 0.5f;
        cxreg[3] = (xlreg[3] + xhreg[3]) * 0.5f;

        cyreg[0] = (ylreg[0] + yhreg[0]) * 0.5f;
        cyreg[1] = (ylreg[1] + yhreg[1]) * 0.5f;
        cyreg[2] = (ylreg[2] + yhreg[2]) * 0.5f;
        cyreg[3] = (ylreg[3] + yhreg[3]) * 0.5f;

        Sample2<FT>(cxreg, cyreg, 0.0f, 0.0f, vreg, 4);
    } // End of binary search.

    if (0 < rn)
        x[m] = cxreg[0];
    if (1 < rn)
        x[m + 1] = cxreg[1];
    if (2 < rn)
        x[m + 2] = cxreg[2];

    if (0 < rn)
        y[m] = cyreg[0];
    if (1 < rn)
        y[m + 1] = cyreg[1];
    if (2 < rn)
        y[m + 2] = cyreg[2];
}

void ConstructQEF2(
    const float* ix0, const float* iy0,
    const float* ix1, const float* iy1,
    const float* nx0, const float* ny0,
    const float* nx1, const float* ny1,
    const int* ens0, const int* ens1, int n,
    float* f, int* h, int* m) {
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
        int c4 = count_horizontal[ens0[i + 1]];

        if (c0 + c1 + c3 < 2){
            h[i] = -1;
            start0 += c0;
            start1 += c2;
            continue;
        }

        h[i] = *m;

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
            const float iy = iy1[start1 + j];
            A[(c0 + j) * 2] = nx;
            A[(c0 + j) * 2 + 1] = ny;
            b[c0 + j] = nx * ix + ny * iy;
            gx += ix / (c0 + c1 + c3);
            gy += iy / (c0 + c1 + c3);
        }

        for (int j = 0; j < c3; ++j) {
            const float nx = nx0[start0 + c0 + c4 + j];
            const float ny = ny0[start0 + c0 + c4 + j];
            const float ix = ix0[start0 + c0 + c4 + j];
            const float iy = iy0[start0 + c0 + c4 + j];
            A[(c0 + c1 + j) * 2] = nx;
            A[(c0 + c1 + j) * 2 + 1] = ny;
            b[c0 + c1 + j] = nx * ix + ny * iy;
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

        memcpy(f + *m * 7, data, sizeof(data));
        *m += 1;

        start0 += c0;
        start1 += c2;
    }
}

void SolveQEF2(const float* f, float* p, int n) {
    const float* s = f;
    float* q = p;
    for (int i = 0; i < n; ++i, s += 7, q += 2) {
        const float ATA[4] = {s[0], s[1], s[1], s[2]};
        float pinv[4];
        PseudoInverse2(ATA, pinv);
        const float ATAg[2] = {
            ATA[0] * s[5] + ATA[1] * s[6],
            ATA[2] * s[5] + ATA[3] * s[6]
        };
        const float d[2] = {
            s[3] - ATAg[0],
            s[4] - ATAg[1]
        };
        const float c[2] = {
            pinv[0] * d[0] + pinv[1] * d[1],
            pinv[2] * d[0] + pinv[3] * d[1]
        };
        q[0] = c[0] + s[5];
        q[1] = c[1] + s[6];
    }
}

inline int NEXTS(int r) {
    return (r + 1) % 3;
}
inline int NNEXTS(int r) {
    return (r + 2) % 3;
}
inline int NEXTR(int r) {
    return (r + 1) % 2;
}
inline int CURR(int r) {
    return r % 2;
}

template <FunctionType FT>
bool GenericDualContour2(
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap, int* pcnt,
    int* edges, int ecap, int* ecnt) {
    int nsamples = n + 1;
    int ncells = n;
    const float xstep = (xt - xs) / (n - 1);
    const float ystep = (yt - ys) / (n - 1);

    // Allocate 3 rows of samples.
    float* const S = (float*)malloc(sizeof(float) * nsamples * 3);
    float* const s[3] = {
        S,
        S + nsamples,
        S + nsamples * 2
    };

    // Initialize two rows of samples.
    Sample2<FT>(
        xs, xt + xstep,
        ys, ys,
        s[0], nsamples);
    Sample2<FT>(
        xs, xt + xstep,
        ys + ystep, ys + ystep,
        s[1], nsamples);

    // Allocate arrays for xlow, ylow, xhigh, yhigh.
    // One cell contains 2 intersected edges at most.
    float* const E = (float*)malloc(sizeof(float) * ncells * 2 * 4);
    float* const xlow = E;
    float* const ylow = E + ncells * 2;
    float* const xhigh = ylow + ncells * 2;
    float* const yhigh = xhigh + ncells * 2;

    // Allocate 4 arrays for ix0, ix1, iy0, iy1.
    // One cell contains 2 intersections at most.
    float* const I = (float*)malloc(sizeof(float) * ncells * 2 * 4);
    float* const ix[2] = {
        I,
        I + ncells * 2 * 2
    };
    float* const iy[2] = {
        I + ncells * 2,
        I + ncells * 2 * 3
    };

    // Allocate 4 arrays for nx0, nx1, ny0, ny1.
    // One cell contains 2 intersections at most.
    float* const N = (float*)malloc(sizeof(float) * ncells * 2 * 4);
    float* const nx[2] = {
        N,
        N + ncells * 2 * 2
    };
    float* const ny[2] = {
        N + ncells * 2,
        N + ncells * 2 * 3
    };

    // Allocate arrays for 2 rows of intersected edges count.
    int* const ENS = (int*)malloc(sizeof(int) * ncells * 2);
    int* const ens[2] = {
        ENS,
        ENS + ncells
    };
    int en[2];

    // Allocate arrays for 2 rows of index of generated mesh vertices.
    int* const H = (int*)malloc(sizeof(int) * ncells * 2);
    int* const h[2] = {
        H,
        H + ncells
    };

    // Allocate arrays for generated mesh vertices.
    float* const V = (float*)malloc(sizeof(float) * 2 * ncells);

    // Allocate buffers for QEF coefficients.
    const int cap = 4 * ncells;
    float* const QEF = (float*)malloc(sizeof(float) * 7 * cap);

    // Initialize the first row of intersections.
    CollectIntersectionEdges2(
        xs, xt + xstep, ys, ys,
        xs, xt + xstep, ys + ystep, ys + ystep,
        s[0], s[1], nsamples,
        xlow, ylow, xhigh, yhigh,
        ens[0], &en[0]);

    SolveIntersection2<FT>(
        xlow, ylow,
        xhigh, yhigh,
        ix[0], iy[0], en[0]);

    SampleGradient2<FT>(
        ix[0], iy[0],
        nx[0], ny[0],
        en[0]);

    float prevy = ys + ystep;
    float y = prevy + ystep;
    // Last count of generated vertices.
    int last = 0;
    *pcnt = 0;
    *ecnt = 0;
    for (int l = 0; l < n - 1; ++l, prevy = y, y += ystep) {
        Sample2<FT>(
            xs, xt + xstep,
            y, y,
            s[NNEXTS(l)], nsamples);

        CollectIntersectionEdges2(
            xs, xt + xstep, prevy, prevy,
            xs, xt + xstep, y, y,
            s[NEXTS(l)], s[NNEXTS(l)], nsamples,
            xlow, ylow, xhigh, yhigh,
            ens[NEXTR(l)], &en[NEXTR(l)]);

        SolveIntersection2<FT>(
            xlow, ylow,
            xhigh, yhigh,
            ix[NEXTR(l)], iy[NEXTR(l)], en[NEXTR(l)]);

        SampleGradient2<FT>(
            ix[NEXTR(l)], iy[NEXTR(l)],
            nx[NEXTR(l)], ny[NEXTR(l)],
            en[NEXTR(l)]);

        int nQEF;
        ConstructQEF2(
            ix[CURR(l)], iy[CURR(l)],
            ix[NEXTR(l)], iy[NEXTR(l)],
            nx[CURR(l)], ny[CURR(l)],
            nx[NEXTR(l)], ny[NEXTR(l)],
            ens[CURR(l)], ens[NEXTR(l)],
            ncells,
            QEF, h[NEXTR(l)], &nQEF);

        if (*pcnt + nQEF > pcap)
            return false;
        SolveQEF2(QEF, p + 2 * (*pcnt), nQEF);
        *pcnt += nQEF;

        // Generate horizontal edges.
        int base1 = *pcnt - nQEF;
        int* const h1 = h[NEXTR(l)];
        for (int c = 1; c < ncells - 1; ++c) {
            if (h1[c] != -1 && h1[c - 1] != -1) {
                if (*ecnt + 1 > ecap)
                    return false;
                edges[*ecnt * 2] = base1 + h1[c];
                edges[*ecnt * 2 + 1] = base1 + h1[c - 1];
                *ecnt += 1;
            }
        }

        // Generate vertical edges.
        int base0 = base1 - last;
        int* const h0 = h[CURR(l)];
        if (l != 0) {
            for (int c = 0; c < ncells - 1; ++c) {
                if (h1[c] != -1 && h0[c] != -1) {
                    if (*ecnt + 1 > ecap)
                        return false;
                    edges[*ecnt * 2] = base1 + h1[c];
                    edges[*ecnt * 2 + 1] = base0 + h0[c];
                    *ecnt += 1;
                }
            }
        }

        last = nQEF;
    }

    free(S);
    free(I);
    free(N);
    free(E);
    free(ENS);
    free(H);
    free(V);
    free(QEF);
    return true;
}

bool DualContour2(
    FunctionType ft,
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap, int* pcnt,
    int* edges, int ecap, int* ecnt) {
    switch(ft) {
        case FT_UNIT_SPHERE:
            return GenericDualContour2<FT_UNIT_SPHERE>(
                xs, xt, ys, yt, n, p, pcap, pcnt, edges, ecap, ecnt);

        case FT_HEART:
            return GenericDualContour2<FT_HEART>(
                xs, xt, ys, yt, n, p, pcap, pcnt, edges, ecap, ecnt);
            
        default:
            assert(false);
            return false;
    }
}
