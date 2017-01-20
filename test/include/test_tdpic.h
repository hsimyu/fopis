#include <tdpic.h>
#include <gtest/gtest.h>

namespace TEST_TDPIC {
    class TDPICTest1 : public ::testing::Test {
        protected:
            TDPICTest1();
            ~TDPICTest1();

            // members
            Environment* env;
            ParticleType* ptype;
            Grid* root_grid;
    };
}
