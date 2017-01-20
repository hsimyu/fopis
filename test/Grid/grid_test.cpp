#include <test_tdpic.h>

using std::cout;
using std::endl;

namespace TEST_TDPIC {
    TEST_F(TDPICTest1, gridGetIDMap)
    {
        // Level 2まで
        root_grid->makeChild(2, 2, 2, 5, 5, 5);
        root_grid->makeChild(4, 8, 8, 9, 9, 9);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);

        std::vector< std::vector<int> > testChild{{0}, {1, 2}, {3}};
        ASSERT_EQ(testChild, root_grid->getIDMapOnRoot());
    }

    TEST_F(TDPICTest1, gridGetChildMap)
    {
        // Level 2まで
        root_grid->makeChild(2, 2, 2, 5, 5, 5);
        root_grid->makeChild(4, 8, 8, 9, 9, 9);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);

        std::map<int, std::vector<int> > testChild;
        testChild[0] = {1, 2};
        testChild[1] = {3};
        testChild[2] = {};
        testChild[3] = {};

        ASSERT_EQ(testChild, root_grid->getChildMapOnRoot());
    }
}
