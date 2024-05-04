/*
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include "utils.h" 


TEST(FindIndex, Test) {
    std::vector<int> testVec = {1, 2, 3, 4, 5};
    EXPECT_EQ(2, clipperplus::find_index(testVec, 3));
    EXPECT_EQ(-1, clipperplus::find_index(testVec, 6));
}

TEST(AdjacencyMatrixToList, SmallGraph) {
    Eigen::MatrixXd adj(3, 3);
    adj << 0, 1, 1,
           1, 0, 1,
           1, 1, 0;
    const int nnodes = 3;
    std::vector<int> ei, ej;

    clipperplus::adjmat_to_adjlist(adj, nnodes, ei, ej);

    // Expected edges: (0,1), (0,2), (1,2)
    ASSERT_EQ(ei[0], 0);
    ASSERT_EQ(ej[0], 1);
    ASSERT_EQ(ei[1], 0);
    ASSERT_EQ(ej[1], 2);
    ASSERT_EQ(ei[2], 1);
    ASSERT_EQ(ej[2], 2);
}

TEST(AdjacencyMatrixToList, InvalidInputs) {
    Eigen::MatrixXd adj(2, 3); // Non-square matrix
    adj << 0, 1, 0,
           1, 0, 1;
    const int nnodes = 2;
    std::vector<int> ei, ej;

    ASSERT_THROW(clipperplus::adjmat_to_adjlist(adj, nnodes, ei, ej), std::invalid_argument);
}

*/
