/**
 * @file clipper.h
 * @brief CLIPPER data association framework
 * @author Parker Lusk <plusk@mit.edu>
 * @date 3 October 2020
 */

#pragma once

#include <tuple>
#include <cmath>

#include <Eigen/Dense>

#include "clipper/invariants/abstract.h"
#include "clipper/invariants/builtins.h"
#include "clipper/types.h"

#include "clipper/dsd.h"
#include "clipper/maxclique.h"

namespace clipper {

/**
 * @brief      CLIPPER parameters
 */
struct Params {
  // gradient descent parameters
  double tol_u = 1e-8; // (1e-8) stop innerloop when change in u < tol
  double tol_F = 1e-9; // (1e-9) stop innerloop when change in F < tol
  int maxiniters = 200; // max num of gradient ascent steps for each d
  int maxoliters = 1000; // max num of outer loop iterations to find d
  // line search parameters
  double beta = 0.5; // backtracking step size reduction, in (0, 1)
  double sigma = 0.01; // threshold of armijo condition
  
  /* maximum number of line search iters per grad step 
    note: do NOT choose large number (>30), or otherwise 
    alpha stepsize will become zero and line search will get stuck */
  int maxlsiters = 20; 
  double minalpha = 1e-9; // minimum value of alpha for line search
  double maxalpha = 2; // maximum value of alpha for line search

  double eps = 1e-9; // (1e-9) numerical threshold around 0 (default 1e-9)

  double affinityeps = 1e-4; ///< sparsity-promoting threshold for affinities

  bool rescale_u0 = true; ///< Rescale u0 using one power iteration. This
                          ///< removes some randomness of the initial guess;
                          ///< i.e., after one step of power method, random
                          ///< u0's look similar.

  // brief Rounding procedure
  enum class Rounding { NONZERO, DSD, DSD_HEU };
  /*  NONZERO - any nonzero elements of u are selected as nodes
      DSD - select the densest edge-weighted subgraph of the
        subgraph induced by NONZERO rounding
      DSD_HEU - A heuristic for selecting the top best nodes
        of the subgraph induced by NONZERO rounding
        DSD_HEU tends to pick smaller subgraphs than DSD,
        sometimes leading to higher precision at the cost
        of lower recall */
  Rounding rounding = Rounding::DSD_HEU;
};

/**
 * @brief      Data associated with a CLIPPER dense clique solution
 */
struct Solution
{
  double t; ///< duration spent solving [s]
  long ifinal; ///< number of outer iterations before convergence
  std::vector<long> nodes; ///< indices of graph vertices in dense clique
  Eigen::VectorXd u0; ///< initial vector used for local solver
  Eigen::VectorXd u; ///< characteristic vector associated with graph
  double score; ///< value of objective function / largest eigenvalue
};

/**
 * @brief      Convenience class to use CLIPPER for data association.
 */
class CLIPPER
{
public:
  CLIPPER(const invariants::PairwiseInvariantPtr& invariant, const Params& params);
  ~CLIPPER() = default;

  /**
 * @brief      Creates an affinity matrix containing consistency scores for
 *             each of the m pairwise associations listed in matrix A.
 *
 * @param[in]  D1           Dataset 1 of n1 d-dim elements (dxn1)
 * @param[in]  D2           Dataset 2 of n2 d-dim elements (dxn2)
 * @param[in]  A            Associations to score (mx2)
 */
  void scorePairwiseConsistency(const invariants::Data& D1,
                                const invariants::Data& D2,
                                const Association& A = Association());

  /**
   * @brief      Solves the MSRC problem using
   * graduated projected gradient ascent
   *
   * @param[in]  u0    Initial condition, if none provided random vec is used
   */
  void solve(const Eigen::VectorXd& u0 = Eigen::VectorXd());

  // find maximal clique in unweighted graph
  void solveBinary(const Eigen::VectorXd& u0 = Eigen::VectorXd());

  /* find maximal clique in unweighted graph
    the old version (identical to matlab) */
  void solveBinary_old(const Eigen::VectorXd& u0 = Eigen::VectorXd());


  /**
   * @brief      Solves the maximum clique problem
   * @param[in]  params Clique solver parameters
   */

  /**
   * @brief      Solves the maximum spectral radius clique problem using a semidefinite relaxation.
   * @param[in]  params  The parameters
   */

  const Solution& getSolution() const { return soln_; }
  Affinity getAffinityMatrix();
  Constraint getConstraintMatrix();

  /**
   * @brief      Skip using scorePairwiseConsistency and directly set the
   *             affinity and constraint matrices. Note that this function
   *             accepts dense matrices. Use the sparse version for better`
   *             performance if you already have sparse matrices available.
   * @param[in]  M     Affinity matrix
   * @param[in]  C     Constraint matrix
   */
  void setMatrixData(const Affinity& M, const Constraint& C);

  void setMatrixData_old(const Affinity& M, const Constraint& C);
  

  /**
   * @brief      Skip using scorePairwiseConsistency and directly set the
   *             affinity and constraint matrices. Note that this function
   *             accepts sparse matrices. These matrices should be upper
   *             triangular and should not have diagonal values set.
   * @param[in]  M     Affinity matrix
   * @param[in]  C     Constraint matrix
   */
  void setSparseMatrixData(const SpAffinity& M, const SpConstraint& C);

  Association getInitialAssociations() const;
  Association getSelectedAssociations() const;

  void setParallelize(bool parallelize) { parallelize_ = parallelize; };

private:
  Params params_;
  invariants::PairwiseInvariantPtr invariant_;

  bool parallelize_ = true; ///< should affinity calculation be parallelized

  Association A_; ///< initial (putative) set of associations

  // \brief Problem data from latest instance of data association
  Solution soln_; ///< solution information from CLIPPER dense clique solver
  SpAffinity M_; ///< affinity matrix (i.e., weighted consistency graph)
  SpConstraint C_; ///< constraint matrix (i.e., prevents forming links)

  Eigen::MatrixXd _M;
  Eigen::MatrixXd _C;

  /**
   * @brief      Identifies a dense clique of an undirected graph G from its
   *             weighted affinity matrix M while satisfying any active
   *             constraints in C (indicated with zeros).
   *
   *             If M is binary and C==M then CLIPPER returns a maximal clique.
   *
   *             This algorithm employs a projected gradient descent method to
   *             solve a symmetric rank-one nonnegative matrix approximation.
   *
   * @param[in]  M        Symmetric, non-negative nxn affinity matrix where
   *                      each element is in [0,1]. Nodes can also be weighted
   *                      between [0,1] (e.g., if there is a prior indication
   *                      that a node belongs to the desired cluster). In the
   *                      case that all nodes are equally likely to be in the
   *                      densest cluster/node weights should not be considered
   *                      set the diagonal of M to identity.
   * @param[in]  C        nxn binary constraint matrix. Active const. are 0.
   */
  void findDenseClique(const Eigen::VectorXd& u0);

  // find maximal clique for a binray adjacency matrix (unweighted graph)
  void findMaximalClique(const Eigen::VectorXd& u0);

  /** find maximal clique for a binary adjacency matrix (unweighted graph)
  * the old version (identical to matlab) */
  void findMaximalClique_old(const Eigen::VectorXd& u0);
};

}