#include <gtest/gtest.h>
#include <dc2-internal.h>
#include <cmath>
using namespace std;

TEST(UnitSphere, Sample2) {
    Function f;
    f.mFunctionType = FT_UNIT_SPHERE;

    float x[3] = {-2.0f, 0.0f, 2.0f};
    float y0[3] = {-2.0f, -2.0f, -2.0f};
    float y1[3] = {0.0f, 0.0f, 0.0f};

    float s[3];

    Sample2<FT_UNIT_SPHERE>(f, x, y0, s, 3);
    EXPECT_FLOAT_EQ(7.0f, s[0]);
    EXPECT_FLOAT_EQ(3.0f, s[1]);
    EXPECT_FLOAT_EQ(7.0f, s[2]);

    Sample2<FT_UNIT_SPHERE>(f, x, y1, s, 3);
    EXPECT_FLOAT_EQ(3.0f, s[0]);
    EXPECT_FLOAT_EQ(-1.0f, s[1]);
    EXPECT_FLOAT_EQ(3.0f, s[2]);    
}

TEST(UnitSphere, SampleGradient2) {
    Function f;
    f.mFunctionType = FT_UNIT_SPHERE;

    float x[3] = {-2.0f, 0.0f, 2.0f};
    float y[3] = {-2.0f, -2.0f, -2.0f};

    float ds0[3], ds1[3];

    SampleGradient2<FT_UNIT_SPHERE>(f, x, y, ds0, ds1, 3);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds0[0]);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds1[0]);

    EXPECT_FLOAT_EQ(0.0f, ds0[1]);
    EXPECT_FLOAT_EQ(-1.0f, ds1[1]);

    EXPECT_FLOAT_EQ(1.0f / sqrt(2.0f), ds0[2]);
    EXPECT_FLOAT_EQ(-1.0f / sqrt(2.0f), ds1[2]);
}

TEST(UnitSphere, CollectIntersectionEdges2) {
    Function f;
    f.mFunctionType = FT_UNIT_SPHERE;

    float x[3] = {-2.0f, 0.0f, 2.0f};
    float y0[3] = {-2.0f, -2.0f, -2.0f};
    float y1[3] = {0.0f, 0.0f, 0.0f};

    float s0[3], s1[3];

    Sample2<FT_UNIT_SPHERE>(f, x, y0, s0, 3);
    Sample2<FT_UNIT_SPHERE>(f, x, y1, s1, 3);

    float xlow[4], ylow[4];
    float xhigh[4], yhigh[4];
    int ens[2];
    int n;

    CollectIntersectionEdges2<FT_UNIT_SPHERE>(
        f,
        x, y0,
        x, y1,
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
    Function f;
    f.mFunctionType = FT_UNIT_SPHERE;

    float x[3] = {-2.0f, 0.0f, 2.0f};
    float y0[3] = {-2.0f, -2.0f, -2.0f};
    float y1[3] = {0.0f, 0.0f, 0.0f};

    float s0[3], s1[3];

    Sample2<FT_UNIT_SPHERE>(f, x, y0, s0, 3);
    Sample2<FT_UNIT_SPHERE>(f, x, y1, s1, 3);

    float xlow[4], ylow[4];
    float xhigh[4], yhigh[4];
    int ens[2];
    int n;

    CollectIntersectionEdges2<FT_UNIT_SPHERE>(
        f,
        x, y0,
        x, y1,
        s0, s1, 3,
        xlow, ylow, xhigh, yhigh,
        ens, &n);

    float ix[4], iy[4];
    SolveIntersection2<FT_UNIT_SPHERE>(
        f,
        xlow, ylow,
        xhigh, yhigh,
        ix, iy, n);

    EXPECT_FLOAT_EQ(0.0f, ix[0]);
    EXPECT_FLOAT_EQ(-1.0f, iy[0]);
}