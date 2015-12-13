#include <gtest/gtest.h>
#include <pinv.h>

#define EXPECT_FLOAT2X2_EQ(a, b)\
    EXPECT_FLOAT_EQ((a)[0], (b)[0]);\
    EXPECT_FLOAT_EQ((a)[1], (b)[1]);\
    EXPECT_FLOAT_EQ((a)[2], (b)[2]);\
    EXPECT_FLOAT_EQ((a)[3], (b)[3])

TEST(Simple, Simple) {
    const float M[] = {1.0f, -2.0f, 3.0f, -1.0f};
    const float E[] = {-0.2f, 0.4f, -0.6f, 0.2f};
    float pinv[4];

    PseudoInverse2(M, pinv);

    EXPECT_FLOAT2X2_EQ(pinv, E);
}
