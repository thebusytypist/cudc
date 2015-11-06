#include <gtest/gtest.h>
#include <dc2-internal.h>
#include <cmath>
using namespace std;

TEST(UnitSphere, Sample2) {
    float s[3];

    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, -2.0f, -2.0f, s, 3);
    EXPECT_FLOAT_EQ(7.0f, s[0]);
    EXPECT_FLOAT_EQ(3.0f, s[1]);
    EXPECT_FLOAT_EQ(7.0f, s[2]);

    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, 0.0f, 0.0f, s, 3);
    EXPECT_FLOAT_EQ(3.0f, s[0]);
    EXPECT_FLOAT_EQ(-1.0f, s[1]);
    EXPECT_FLOAT_EQ(3.0f, s[2]);    
}

TEST(UnitSphere, SampleGradient2) {
    float x[3] = {-2.0f, 0.0f, 2.0f};
    float y[3] = {-2.0f, -2.0f, -2.0f};

    float ds0[3], ds1[3];

    SampleGradient2<FT_UNIT_SPHERE>(x, y, ds0, ds1, 3);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds0[0]);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds1[0]);

    EXPECT_FLOAT_EQ(0.0f, ds0[1]);
    EXPECT_FLOAT_EQ(-1.0f, ds1[1]);

    EXPECT_FLOAT_EQ(1.0f / sqrt(2.0f), ds0[2]);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds1[2]);
}

TEST(UnitSphere, CollectIntersectionEdges2) {
    float s0[3], s1[3];

    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, -2.0f, -2.0f, s0, 3);
    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, 0.0f, 0.0f, s1, 3);

    float xlow[4], ylow[4];
    float xhigh[4], yhigh[4];
    int ens[2];
    int n;

    CollectIntersectionEdges2(
        -2.0f, 2.0f, -2.0f, -2.0f,
        -2.0f, 2.0f, 0.0f, 0.0f,
        s0, s1, 3,
        xlow, ylow, xhigh, yhigh,
        ens, &n);

    EXPECT_EQ(1, n);
    EXPECT_EQ(0, ens[0]);
    EXPECT_EQ(2, ens[1]);

    EXPECT_FLOAT_EQ(0.0f, xlow[0]);
    EXPECT_FLOAT_EQ(0.0f, ylow[0]);
    EXPECT_FLOAT_EQ(0.0f, xhigh[0]);
    EXPECT_FLOAT_EQ(-2.0f, yhigh[0]);
}

TEST(UnitSphere, SolveIntersection2) {
    float s0[3], s1[3];

    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, -2.0f, -2.0f, s0, 3);
    Sample2<FT_UNIT_SPHERE>(-2.0f, 2.0f, 0.0f, 0.0f, s1, 3);

    float xlow[4], ylow[4];
    float xhigh[4], yhigh[4];
    int ens[2];
    int n;

    CollectIntersectionEdges2(
        -2.0f, 2.0f, -2.0f, -2.0f,
        -2.0f, 2.0f, 0.0f, 0.0f,
        s0, s1, 3,
        xlow, ylow, xhigh, yhigh,
        ens, &n);

    float ix[4], iy[4];
    SolveIntersection2<FT_UNIT_SPHERE>(
        xlow, ylow,
        xhigh, yhigh,
        ix, iy, n);

    EXPECT_FLOAT_EQ(0.0f, ix[0]);
    EXPECT_FLOAT_EQ(-1.0f, iy[0]);
}

TEST(UnitSphere, QEF2) {
    float s0[4], s1[4], s2[4], s3[4];

    Sample2<FT_UNIT_SPHERE>(-2.0f, 4.0f, -2.0f, -2.0f, s0, 4);
    Sample2<FT_UNIT_SPHERE>(-2.0f, 4.0f, 0.0f, 0.0f, s1, 4);
    Sample2<FT_UNIT_SPHERE>(-2.0f, 4.0f, 2.0f, 2.0f, s2, 4);
    Sample2<FT_UNIT_SPHERE>(-2.0f, 4.0f, 4.0f, 4.0f, s3, 4);

    float xlow[6], ylow[6];
    float xhigh[6], yhigh[6];
    int ens0[3], ens1[3];
    int n0, n1;
    float ix0[6], iy0[6], ix1[6], iy1[6];
    float nx0[6], ny0[6];
    float nx1[6], ny1[6];

    CollectIntersectionEdges2(
        -2.0f, 4.0f, -2.0f, -2.0f,
        -2.0f, 4.0f, 0.0f, 0.0f,
        s0, s1, 4,
        xlow, ylow, xhigh, yhigh,
        ens0, &n0);
    SolveIntersection2<FT_UNIT_SPHERE>(
        xlow, ylow,
        xhigh, yhigh,
        ix0, iy0,
        n0);
    SampleGradient2<FT_UNIT_SPHERE>(
        ix0, iy0,
        nx0, ny0,
        n0);

    CollectIntersectionEdges2(
        -2.0f, 4.0f, 0.0f, 0.0f,
        -2.0f, 4.0f, 2.0f, 2.0f,
        s1, s2, 4,
        xlow, ylow, xhigh, yhigh,
        ens1, &n1);
    SolveIntersection2<FT_UNIT_SPHERE>(
        xlow, ylow,
        xhigh, yhigh,
        ix1, iy1,
        n1);
    SampleGradient2<FT_UNIT_SPHERE>(
        ix1, iy1,
        nx1, ny1,
        n1);
    
    float qef[7 * 2];
    int h[2];
    int m = 0;
    ConstructQEF2(
        ix0, iy0,
        ix1, iy1,
        nx0, ny0,
        nx1, ny1,
        ens0, ens1, 3,
        qef, h, &m);
    EXPECT_EQ(0, h[0]);
    EXPECT_EQ(1, h[1]);

    float p[2 * 2];
    SolveQEF2(qef, p, m);

    EXPECT_EQ(2, m);
    EXPECT_FLOAT_EQ(-1.0f, p[0]);
    EXPECT_FLOAT_EQ(-1.0f, p[1]);
    EXPECT_FLOAT_EQ(1.0f, p[2]);
    EXPECT_FLOAT_EQ(-1.0f, p[3]);

    // Check the second row.
    CollectIntersectionEdges2(
        -2.0f, 4.0f, 2.0f, 2.0f,
        -2.0f, 4.0f, 4.0f, 4.0f,
        s2, s3, 4,
        xlow, ylow, xhigh, yhigh,
        ens0, &n0);
    SolveIntersection2<FT_UNIT_SPHERE>(
        xlow, ylow,
        xhigh, yhigh,
        ix0, iy0,
        n0);
    SampleGradient2<FT_UNIT_SPHERE>(
        ix0, iy0,
        nx0, ny0,
        n0);

    ConstructQEF2(
        ix1, iy1,
        ix0, iy0,
        nx1, ny1,
        nx0, ny0,
        ens1, ens0, 3,
        qef, h, &m);
    EXPECT_EQ(0, h[0]);
    EXPECT_EQ(1, h[1]);

    SolveQEF2(qef, p, m);

    EXPECT_EQ(2, m);
    EXPECT_FLOAT_EQ(-1.0f, p[0]);
    EXPECT_FLOAT_EQ(1.0f, p[1]);
    EXPECT_FLOAT_EQ(1.0f, p[2]);
    EXPECT_FLOAT_EQ(1.0f, p[3]);
}

TEST(UnitSphere, DualContour2) {
    const int pcap = 4;
    int pcnt = 0;
    float p[2 * pcap];

    const int ecap = 4;
    int ecnt = 0;
    int edges[2 * ecap];

    bool r = DualContour2(FT_UNIT_SPHERE,
        -2.0f, 2.0f,
        -2.0f, 2.0f, 3,
        p, pcap, &pcnt,
        edges, ecap, &ecnt);

    EXPECT_TRUE(r);

    // Check vertices.
    EXPECT_EQ(4, pcnt);
    float ep[8] = {
        -1.0f, -1.0f,
        1.0f, -1.0f,
        -1.0f, 1.0f,
        1.0f, 1.0f
    };
    for (int i = 0; i < 8; ++i) {
        EXPECT_FLOAT_EQ(ep[i], p[i]);
    }

    // Check edges.
    EXPECT_EQ(4, ecnt);
    int ee[8] = {
        1, 0,
        3, 2,
        2, 0,
        3, 1
    };
    for (int i = 0; i < 8; ++i) {
        EXPECT_EQ(ee[i], edges[i]);
    }
}
