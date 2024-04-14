/*
c++ tests for clipper+ library
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "clipperplus/utils.h"

// for computating core numbers and PMC-heuristic clique 
#include "clipperplus/clipperplus_clique.h"

TEST(CLIPPERPLUS, clique1) {
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

    long clique_size = 0;
    std::vector<int> clique;
    Certificate certificate = NONE;
    clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);

    EXPECT_EQ(clique_size, 6);
    std::vector<int> clique_expected = {1, 3, 4, 6, 7, 8};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, NONE);
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

    long clique_size = 0;
    std::vector<int> clique;
    Certificate certificate = NONE;
    clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);

    EXPECT_EQ(clique_size, 7);
    std::vector<int> clique_expected = {6, 0, 2, 3, 5, 8, 9};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, CHROMATIC);
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
    long clique_size = 0;
    std::vector<int> clique;
    Certificate certificate = NONE;
    clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);

    EXPECT_EQ(clique_size, 8);
    std::vector<int> clique_expected = {4, 10, 13, 14, 15, 16, 17, 18};
    EXPECT_EQ(clique, clique_expected); 
    EXPECT_EQ(certificate, NONE);
};
