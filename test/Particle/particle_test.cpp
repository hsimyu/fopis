#include <tdpic.h>
#include <gtest/gtest.h>

namespace {
    int add(int x, int y)
    {
        return x + y;
    }

    TEST(AddTest, Test1)
    {
        ASSERT_EQ(2, add(1, 1));
        ASSERT_EQ(10, add(-1, 11));
        ASSERT_EQ(10, add(0, 10));
    }
}
