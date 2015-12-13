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

TEST(Simple, FloatAccuracy1) {
    const float M[] = {0.621373f, 0.925549f, 0.925549f, 1.378627f};
    const float E[] = {
        2142488.3841422f, -1438371.64182511f,
        -1438371.64182511f, 965659.62709245f};
    float pinv[4];
    PseudoInverse2(M, pinv);
    
    EXPECT_FLOAT2X2_EQ(pinv, E);
}

TEST(Simple, FloatAccuracy2) {
    const float M[] = {1.560416f, 0.828206f, 0.828206f, 0.439584f};
    const float E[] = {
        50361.87169762f, -94885.17396264f,
        -94885.17396264f, 178772.36293156f};
    float pinv[4];
    PseudoInverse2(M, pinv);

    EXPECT_FLOAT2X2_EQ(pinv, E);
}
