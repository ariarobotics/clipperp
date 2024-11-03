/**
 * @file wrappers.h
 * @brief Wrapper for Clipperplus' C++ functions with parameters passed by reference
 */

#pragma once

#include <pybind11/pybind11.h>
#include "clipperplus/clipperplus_clique.h"
#include "clipperplus/clique_optimization.h"

class Wrapper {
  public:
    static std::tuple<long, std::vector<int>, int> clipperplus_clique_wrapper(const Eigen::MatrixXd& adj){
      auto [clique, certificate] = clipperplus::find_clique(adj);

      return std::make_tuple((long)clique.size(), clique, (int)certificate);
    }

    static std::vector<int> find_heuristic_clique_wrapper(
      const Eigen::MatrixXd& adj, 
      std::vector<int>& clique
    ){
      clique = clipperplus::find_heuristic_clique(adj);
      return clique;
    }
    
    static std::tuple<int, unsigned long, std::vector<long>> clique_optimization_wrapper(
      const Eigen::MatrixXd& M, 
      const Eigen::VectorXd& u0
    ){
      std::vector<long> clique = clipperplus::clique_optimization(M, u0, clipperplus::Params());
      return std::make_tuple(1, clique.size(), clique);
    }


    static std::vector<int> tcore_heuristic(const Eigen::MatrixXd& adj) {
      std::vector<int> triangle_core = clipperplus::Graph(adj).get_triangle_core();

      int max_tcore = *std::max_element(triangle_core.begin(), triangle_core.end());
      std::vector<int> heuristic_clique;
      for (int i = 0; i < triangle_core.size(); ++i) {
        if (triangle_core[i] == max_tcore) {
          heuristic_clique.push_back(i);
        }
      }

      return heuristic_clique;
    }

};

