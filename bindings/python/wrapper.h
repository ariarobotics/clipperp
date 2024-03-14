/**
 * @file wrappers.h
 * @brief Wrapper for Clipperplus' C++ functions with parameters passed by reference
 */

#pragma once

#include <pybind11/pybind11.h>
#include "clipperplus/clipperplus_clique.h"
#include "clipperplus/clique_corenumber.h"
#include "clipperplus/clique_optimization.h"

class Wrapper {
  public:
    static std::tuple<long, std::vector<int>, int> clipperplus_clique_wrapper(const Eigen::MatrixXd& adj){
      long clique_size = 0;
      std::vector<int> clique;
      int certificate = 0;
      clipperplus::clipperplus_clique(adj, clique_size, clique, certificate);
      return std::make_tuple(clique_size, clique, certificate);
    }

    static std::tuple<int, std::vector<int>> find_heuristic_clique_wrapper(const Eigen::MatrixXd& adj, 
                                                                                const std::vector<long>& core_numbers, 
                                                                                std::vector<int>& clique){
      int output = clipperplus::find_heuristic_clique(adj, core_numbers, clique);
      return std::make_tuple(output, clique);
    }

    static std::tuple<unsigned long, std::vector<int>, std::vector<long>, long, std::vector<int>, int> clique_corenumber_wrapper(const Eigen::MatrixXd& adj,
                                                                                  std::vector<int>& clique,
                                                                                  std::vector<long>& core_numbers,
                                                                                  std::vector<int>& node_colors){
      long core_bound = 0;
      int chromatic_bound = 0;
      unsigned long output = clipperplus::clique_corenumber(adj, clique, core_numbers, core_bound, node_colors, chromatic_bound);
      return std::make_tuple(output, clique, core_numbers, core_bound, node_colors, chromatic_bound);
    }
    
    static std::tuple<int, unsigned long, std::vector<long>> clique_optimization(const Eigen::MatrixXd& M, 
                                                                const Eigen::VectorXd& u0){
      unsigned long clique_size = 0;
      std::vector<long> clique;
      int output = clipperplus::clique_optimization(M, u0, clique_size, clique);
      return std::make_tuple(output, clique_size, clique);
    }

};

