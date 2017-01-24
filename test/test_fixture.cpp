#include <test_tdpic.h>

using std::cout;
using std::endl;

namespace TEST_TDPIC {
    TDPICTest1::TDPICTest1() {
        Initializer::initTDPIC(root_grid);
    }

    // テスト毎に実行される
    TDPICTest1::~TDPICTest1() {
        Grid::resetNextID();
    }
}
