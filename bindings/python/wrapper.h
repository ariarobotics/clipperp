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
    static std::tuple<long, std::vector<long>, int> clipperplus_clique_wrapper(const Eigen::MatrixXd& adj){
      auto [clique, certificate] = clipperplus::find_clique(adj);

      return std::make_tuple((long)clique.size(), clique, (int)certificate);
    }

    static std::vector<float> kcore_wrapper(const Eigen::MatrixXd& adj){
      clipperplus::Graph graph(adj);
      return graph.get_core_numbers();
    }

    static std::vector<long> find_heuristic_clique_wrapper(
      const Eigen::MatrixXd& adj, 
      std::vector<long>& clique
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

};

