/*
c++ tests for clipper+ library
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>

// for computating core numbers and PMC-heuristic clique 
#include "clipperplus/clipperplus_clique.h"


TEST(CLIPPERPLUS, clique1) 
{
    std::cout << "adjacency matrix 1:\n" << std::endl;

    Eigen::MatrixXd adj(10,10); // graph affinity matrix
    adj << 0, 0, 1, 1, 1, 1, 1, 0, 1, 0,
           0, 0, 1, 1, 1, 0, 1, 1, 1, 1,
           1, 1, 0, 1, 0, 1, 1, 1, 0, 1,
           1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 0, 0, 1, 1, 1, 1,
           1, 0, 1, 1, 0, 0, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 0, 1, 1, 0,
           0, 1, 1, 1, 1, 1, 1, 0, 1, 1,
           1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
           0, 1, 1, 1, 1, 1, 0, 1, 1, 0;

    std::cout << adj << "\n" << std::endl;

    auto [clique, certificate] = clipperplus::find_clique(adj);
    std::sort(clique.begin(), clique.end());

    EXPECT_EQ(clique.size(), 6);
    decltype(clique) clique_expected = {1, 3, 4, 6, 7, 8};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, clipperplus::CERTIFICATE::NONE);

    // in this case we have two equal maximum cliques and it fails!!
    // one should update the test-case
};

TEST(CLIPPERPLUS, clique2) {
    std::cout << "adjacency matrix 2: \n" << std::endl;

    Eigen::MatrixXd adj(10,10); // graph affinity matrix
    adj << 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
           0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
           1, 1, 1, 1, 0, 1, 0, 1, 1, 0,
           1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
           1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
           1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
           1, 1, 1, 1, 0, 1, 1, 1, 1, 0;   

    std::cout << adj << "\n" << std::endl;

    auto [clique, certificate] = clipperplus::find_clique(adj);
    std::sort(clique.begin(), clique.end());

    EXPECT_EQ(clique.size(), 7);
    decltype(clique) clique_expected = {0, 2, 3, 5, 6, 8, 9};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, clipperplus::CERTIFICATE::CHROMATIC_BOUND);
};

TEST(CLIPPERPLUS, clique3) {
    std::cout << "adjacency matrix 3: \n" << std::endl;

    Eigen::MatrixXd adj(20,20); // graph affinity matrix
    adj << 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
           0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1,
           1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0,
           0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1,
           0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
           0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0,
           0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
           1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1,
           0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
           1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
           1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
           1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0,
           1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
           0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
           1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0;

    std::cout << adj << "\n" << std::endl;  

    auto [clique, certificate] = clipperplus::find_clique(adj);
    std::sort(clique.begin(), clique.end());

    EXPECT_EQ(clique.size(), 8);
    decltype(clique) clique_expected = {4, 10, 13, 14, 15, 16, 17, 18};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, clipperplus::CERTIFICATE::NONE);
};


TEST(CLIPPERPLUS, clique4) {
    std::cout << "adjacency matrix 1:\n" << std::endl;

    Eigen::MatrixXd adj(5, 5); // graph affinity matrix
    adj << 0, 1, 1, 0, 0,
           1, 0, 1, 0, 0,
           1, 1, 0, 0, 0,
           0, 0, 0, 0, 1,
           0, 0, 0, 1, 0;

    std::cout << adj << "\n" << std::endl;

    auto [clique, certificate] = clipperplus::find_clique(adj);
    std::sort(clique.begin(), clique.end());

    EXPECT_EQ(clique.size(), 3);
    decltype(clique) clique_expected = {0, 1, 2};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, clipperplus::CERTIFICATE::HEURISTIC);
};
