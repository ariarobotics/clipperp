#include <gtest/gtest.h>
#include <vector>
#include "clique_corenumber.h"


TEST(CliqueCoreNumberTest, CompleteGraph) {
    Eigen::MatrixXd adj(4, 4);
    adj << 0, 1, 1, 1,
           1, 0, 1, 1,
           1, 1, 0, 1,
           1, 1, 1, 0;
    std::vector<int> clique;
    std::vector<long> core_numbers(4);
    long core_bound;
    std::vector<int> node_colors(4);
    int chromatic_bound;

    int clique_size = clipperplus::clique_corenumber(adj, clique, core_numbers, core_bound, node_colors, chromatic_bound);

    EXPECT_EQ(clique_size, 4);
    EXPECT_EQ(core_bound, 4);
    for (const long core_number : core_numbers) {
        EXPECT_EQ(core_number, 3); 
    }
    EXPECT_EQ(chromatic_bound, 4);
}


TEST(CliqueCoreNumberTest, NormalGraph) {
    Eigen::MatrixXd adj(4, 4);
    adj << 0, 1, 1, 1,
           1, 0, 0, 1,
           1, 0, 0, 1,
           1, 1, 1, 0;
    std::vector<int> clique;
    std::vector<long> core_numbers(4);
    long core_bound;
    std::vector<int> node_colors(4);
    int chromatic_bound;

    int clique_size = clipperplus::clique_corenumber(adj, clique, core_numbers, core_bound, node_colors, chromatic_bound);

    EXPECT_EQ(clique_size, 3);

    EXPECT_EQ(core_bound, 3);

    EXPECT_EQ(core_numbers[0], 2); 
    EXPECT_EQ(core_numbers[1], 2); 
    EXPECT_EQ(core_numbers[2], 2); 
    EXPECT_EQ(core_numbers[3], 2); 

    EXPECT_EQ(chromatic_bound, 3);
}

