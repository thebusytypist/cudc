#include <gtest/gtest.h>
#include "cudc2.h"

TEST(CUDC2, Simple) {
    bool r = StreamingCompact();
    EXPECT_TRUE(r);
}
