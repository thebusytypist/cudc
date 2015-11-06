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
void SampleGradient2<FT_UNIT_SPHERE>(
    const float* x, const float* y,
    float* ds0, float* ds1,
    int n) {
    const float dh = 1e-5f;

    for (int i = 0; i < n; ++i) {
        ds0[i] = 2 * x[i];
        ds1[i] = 2 * y[i];
        float l = sqrt(ds0[i] * ds0[i] + ds1[i] * ds1[i]);
        ds0[i] /= l;
        ds1[i] /= l;
    }
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

template <>
void SolveIntersection2<FT_UNIT_SPHERE>(
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n) {
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

bool ConstructQEF2(
    const float* ix0, const float* iy0,
    const float* ix1, const float* iy1,
    const float* nx0, const float* ny0,
    const float* nx1, const float* ny1,
    const int* ens0, const int* ens1, int n,
    float* f, int* h, int* m, int capacity) {
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
        
        if (*m >= capacity)
            return false;

        memcpy(f + i * 7, data, sizeof(data));
        *m += 1;

        start0 += c0;
        start1 += c2;
    }

    return true;
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
    const float ystep = (xt - xs) / (n - 1);

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
    for (int l = 0; l < n - 1; ++l, y += ystep) {
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
        bool result = ConstructQEF2(
            ix[CURR(l)], iy[CURR(l)],
            ix[NEXTR(l)], iy[NEXTR(l)],
            nx[CURR(l)], ny[CURR(l)],
            nx[NEXTR(l)], ny[NEXTR(l)],
            ens[CURR(l)], ens[NEXTR(l)],
            ncells,
            QEF, h[NEXTR(l)], &nQEF, cap);
        if (!result)
            return false;

        if (pcap < nQEF)
            return false;
        SolveQEF2(QEF, p + 2 * (*pcnt), nQEF);
        *pcnt += nQEF;
        pcap -= nQEF;

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
    int* edge, int ecap, int* ecnt) {
    switch(ft) {
        case FT_UNIT_SPHERE:
            return GenericDualContour2<FT_UNIT_SPHERE>(
                xs, xt, ys, yt, n, p, pcap, pcnt, edge, ecap, ecnt);
        default:
            assert(false);
            return false;
    }
}
