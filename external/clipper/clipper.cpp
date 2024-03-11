/**
 * @file clipper.cpp
 * @brief CLIPPER data association framework
 * @author Parker Lusk <plusk@mit.edu>
 * @date 19 March 2022
 */

#include <iostream>

#include "clipper/clipper.h"
#include "clipper/utils.h"

namespace clipper {

CLIPPER::CLIPPER(const invariants::PairwiseInvariantPtr& invariant, const Params& params)
: params_(params), invariant_(invariant)
{}

// ----------------------------------------------------------------------------

void CLIPPER::scorePairwiseConsistency(const invariants::Data& D1,
                              const invariants::Data& D2, const Association& A)
{
  if (A.size() == 0) A_ = utils::createAllToAll(D1.cols(), D2.cols());
  else A_ = A;

  const size_t m = A_.rows();

  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m, m);

#pragma omp parallel for shared(A_, D1, D2, M_, C_) if(parallelize_)
  for (size_t k=0; k<m*(m-1)/2; ++k) {
    size_t i;
    size_t j; 
    std::tie(i, j) = utils::k2ij(k, m);

    if (A_(i,0) == A_(j,0) || A_(i,1) == A_(j,1)) {
      // violates distinctness constraint
      continue;
    }

    //
    // Evaluate the consistency of geometric invariants associated with ei, ej
    //

    // points to extract invariant from in D1
    const auto& d1i = D1.col(A_(i,0));
    const auto& d1j = D1.col(A_(j,0));

    // points to extract invariant from in D2
    const auto& d2i = D2.col(A_(i,1));
    const auto& d2j = D2.col(A_(j,1));

    const double scr = (*invariant_)(d1i, d1j, d2i, d2j);
    if (scr > params_.affinityeps) { // does not violate inconsistency constraint
      M(i,j) = scr;
    }
  }

  // Identity on diagonal is taken care of implicitly in findDenseClique()
  // M += Eigen::MatrixXd::Identity(m, m);

  M_ = M.sparseView();

  C_ = M_;
  C_.coeffs() = 1;
}

// ----------------------------------------------------------------------------

void CLIPPER::solve(const Eigen::VectorXd& _u0)
{
  Eigen::VectorXd u0;
  if (_u0.size() == 0) {
    u0 = utils::randvec(M_.cols());
  } else {
    u0 = _u0;
  }
  findDenseClique(u0);
}


// ----------------------------------------------------------------------------

void CLIPPER::solveBinary(const Eigen::VectorXd& _u0)
{
  Eigen::VectorXd u0;
  if (_u0.size() == 0) {
    u0 = utils::randvec(M_.cols());
  } else {
    u0 = _u0;
  }
  findMaximalClique(u0);
}

// ----------------------------------------------------------------------------

void CLIPPER::solveBinary_old(const Eigen::VectorXd& _u0)
{
  Eigen::VectorXd u0;
  if (_u0.size() == 0) {
    u0 = utils::randvec(M_.cols());
  } else {
    u0 = _u0;
  }
  findMaximalClique_old(u0);
}

// ----------------------------------------------------------------------------

Association CLIPPER::getInitialAssociations() const
{
  return A_;
}

// ----------------------------------------------------------------------------

Association CLIPPER::getSelectedAssociations() const
{
  return utils::selectInlierAssociations(soln_, A_);
}

// ----------------------------------------------------------------------------

Affinity CLIPPER::getAffinityMatrix()
{
  Affinity M = SpAffinity(M_.selfadjointView<Eigen::Upper>())
                + Affinity::Identity(M_.rows(), M_.cols());
  return M;
}

// ----------------------------------------------------------------------------

Constraint CLIPPER::getConstraintMatrix()
{
  Constraint C = SpConstraint(C_.selfadjointView<Eigen::Upper>())
                  + Constraint::Identity(C_.rows(), C_.cols());
  return C;
}

// ----------------------------------------------------------------------------

void CLIPPER::setMatrixData(const Affinity& M, const Constraint& C)
{
  Eigen::MatrixXd MM = M.triangularView<Eigen::Upper>();
  MM.diagonal().setZero();
  M_ = MM.sparseView();

  Eigen::MatrixXd CC = C.triangularView<Eigen::Upper>();
  CC.diagonal().setZero();
  C_ = CC.sparseView();
}

// ----------------------------------------------------------------------------

void CLIPPER::setMatrixData_old(const Affinity& M, const Constraint& C)
{
  _M = M;
  _C = C;
}

// ----------------------------------------------------------------------------

void CLIPPER::setSparseMatrixData(const SpAffinity& M, const SpConstraint& C)
{
  M_ = M;
  C_ = C;
}

// ----------------------------------------------------------------------------
// Private Methods
// ----------------------------------------------------------------------------

void CLIPPER::findDenseClique(const Eigen::VectorXd& u0)
{
  const auto t_start = std::chrono::high_resolution_clock::now();

  //
  // Initialization
  //

  const size_t n = M_.cols();
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);

  // initialize memory
  Eigen::VectorXd gradF(n);
  Eigen::VectorXd gradFnew(n);
  Eigen::VectorXd u(n);
  Eigen::VectorXd unew(n);
  Eigen::VectorXd Mu(n);
  Eigen::VectorXd num(n);
  Eigen::VectorXd den(n);

  // one step of power method to have a good scaling of u
  if (params_.rescale_u0) {
    u = M_.selfadjointView<Eigen::Upper>() * u0 + u0;
  } else {
    u = u0;
  }
  u /= u.norm();

  // initial value of d
  double d = 0; // zero if there are no active constraints
  Eigen::VectorXd Cbu = ones * u.sum() - C_.selfadjointView<Eigen::Upper>() * u - u;
  const Eigen::VectorXi idxD = ((Cbu.array() > params_.eps) && (u.array() > params_.eps)).cast<int>();
  if (idxD.sum() > 0) {
    Mu = M_.selfadjointView<Eigen::Upper>() * u + u;
    num = utils::selectFromIndicator(Mu, idxD);
    den = utils::selectFromIndicator(Cbu, idxD);
    d = (num.array() / den.array()).minCoeff();    
  }
  #ifdef DEBUG_OPTIM
      std::cout << "clipper: initial d: " << d << std::endl;
  #endif

  #ifdef DEBUG_TIMING
      const auto t2 = std::chrono::high_resolution_clock::now(); // timer
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t_start);
      auto elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: u0, d set up: " << elaps << std::endl;
  #endif

  // Orthogonal projected gradient ascent 

  double F = 0; // objective value

  // iteration counters

  size_t i;
  size_t j; 
  for (i=0; i<params_.maxoliters; ++i) {
    gradF = (1 + d) * u - d * ones * u.sum() + M_.selfadjointView<Eigen::Upper>() * u + C_.selfadjointView<Eigen::Upper>() * u * d;
    F = u.dot(gradF); // current objective value

    // Orthogonal projected gradient ascent

    for (j=0; j<params_.maxiniters; ++j) {
      double alpha = 1;

      // Backtracking line search on gradient ascent

      double Fnew = 0;
      double deltaF = 0;
      for (size_t k=0; k<params_.maxlsiters; ++k) {
        unew = u + alpha * gradF;                     // gradient step
        unew = unew.cwiseMax(0);                      // project onto positive orthant
        unew.normalize();                             // project onto S^n
        gradFnew = (1 + d) * unew // because M/C is missing identity on diagonal
                    - d * ones * unew.sum()
                    + M_.selfadjointView<Eigen::Upper>() * unew
                    + C_.selfadjointView<Eigen::Upper>() * unew * d;
        Fnew = unew.dot(gradFnew);                    // new objective value after step

        deltaF = Fnew - F;                            // change in objective value

        if (deltaF < -params_.eps) {

          // objective value decreased---we need to backtrack, so reduce step size

          alpha = alpha * params_.beta;
        } else {
          break; // obj value increased, stop line search
        }

        #ifdef DEBUG_OPTIM
          if (k==params_.maxlsiters-1) {
            std::cout << "clipper: reached max line search itr." << std::endl;
          }
        #endif
      }
      const double deltau = (unew - u).norm();

      // update values
      F = Fnew;
      u = unew;
      gradF = gradFnew;

      // check if desired accuracy has been reached by gradient ascent 
      if (deltau < params_.tol_u || std::abs(deltaF) < params_.tol_F) break;      
    }

    // Increase d

    Cbu = ones * u.sum() - C_.selfadjointView<Eigen::Upper>() * u - u;
    const Eigen::VectorXi idxD = ((Cbu.array() > params_.eps) && (u.array() > params_.eps)).cast<int>();
    if (idxD.sum() > 0) {
      Mu = M_.selfadjointView<Eigen::Upper>() * u + u;
      num = utils::selectFromIndicator(Mu, idxD);
      den = utils::selectFromIndicator(Cbu, idxD);
      const double deltad = (num.array() / den.array()).abs().minCoeff();

      d += deltad;

      #ifdef DEBUG_OPTIM
        std::cout << "clipper: u: ";
        for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";}
        std::cout << std::endl;

        std::cout << "clipper: idxD: ";
        for (int ii=0; ii<idxD.size(); ++ii) {std::cout << idxD(ii) << " ";}
        std::cout << std::endl;
        
        std::cout << "clipper: Mu: "; 
        for (int ii=0; ii<Mu.size(); ++ii) {std::cout << Mu(ii) << " ";}
        std::cout << std::endl;

        std::cout << "clipper: Cbu: ";
        for (int ii=0; ii<Cbu.size(); ++ii) {std::cout << Cbu(ii) << " ";}
        std::cout << std::endl;

        std::cout << "clipper: num: ";
        for (int ii=0; ii<num.size(); ++ii) {std::cout << num(ii) << " ";}
        std::cout << std::endl;

        std::cout << "clipper: den: ";
        for (int ii=0; ii<den.size(); ++ii) {std::cout << den(ii) << " ";}
        std::cout << std::endl;
      #endif

    } else {
      break;
    }

    #ifdef DEBUG_OPTIM
      std::cout << "clipper: j loop: " << j << std::endl;
      std::cout << "clipper: increased d: " << d << std::endl;
    #endif

  }

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: i loop: " << i << std::endl;
  #endif

  #ifdef DEBUG_TIMING
      const auto t3 = std::chrono::high_resolution_clock::now(); // timer
      dur = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
      elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: optimization: " << elaps << std::endl;
  #endif

  // Generate output

  // node indices of rounded u vector
  std::vector<long> nodes;

  #ifdef DEBUG_OPTIM
      std::cout << "clipper: u final:\n" << u << std::endl;
  #endif

  if (params_.rounding == Params::Rounding::NONZERO) {

    nodes = utils::findIndicesWhereAboveThreshold(u, 0.0);    

  } else if (params_.rounding == Params::Rounding::DSD) {

    // subgraph induced by non-zero elements of u
    const std::vector<long> S = utils::findIndicesWhereAboveThreshold(u, 0.0);

    // TODO(plusk): make this faster by leveraging matrix sparsity
    nodes = dsd::solve(M_, S);

  } else if (params_.rounding == Params::Rounding::DSD_HEU) {

    // estimate cluster size using largest eigenvalue
    auto omega = static_cast<int>(std::round(F));

    // extract indices of nodes in identified dense cluster
    nodes = utils::findIndicesOfkLargest(u, omega);

  }
   

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: nodes above thresh: ";
    for (long i : nodes) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  #endif

  #ifdef DEBUG_TIMING
      const auto t4 = std::chrono::high_resolution_clock::now(); // timer
      dur = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
      elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: rounding output: " << elaps << std::endl;
  #endif

  const auto t_end = std::chrono::high_resolution_clock::now();
  const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
  const double elapsed = static_cast<double>(duration.count()) / 1e9;

  // set solution
  soln_.t = elapsed;
  soln_.ifinal = i;
  std::swap(soln_.nodes, nodes);
  soln_.u0 = u0;
  soln_.u.swap(u);
  soln_.score = F;
}



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



// find maximal clique for a binray adjacency matrix (i.e., unweighted graph)
void CLIPPER::findMaximalClique(const Eigen::VectorXd& u0)
{
  #ifdef DEBUG_OPTIM
    std::cout << "running clipper binary version..." << std::endl;
  #endif

  const auto t_start = std::chrono::high_resolution_clock::now();
  
  // Initialization
  const size_t n = M_.cols();
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
  // initialize memory
  Eigen::VectorXd gradF(n);
  Eigen::VectorXd gradFnew(n);
  Eigen::VectorXd u(n);
  Eigen::VectorXd unew(n);
  Eigen::VectorXd uold(n);
  Eigen::VectorXd Mu(n);
  Eigen::VectorXd num(n);
  Eigen::VectorXd den(n);

  // one step of power method to have a good scaling of u
  if (params_.rescale_u0) {
    u = M_.selfadjointView<Eigen::Upper>() * u0 + u0;
  } else {
    u = u0;
  }
  u /= u.norm();

  // initial value of d
  double d = 0; // zero if there are no active constraints
  Eigen::VectorXd Cbu = ones * u.sum() - C_.selfadjointView<Eigen::Upper>() * u - u;
  const Eigen::VectorXi idxD = ((Cbu.array() > params_.eps) && (u.array() > params_.eps)).cast<int>();
  if (idxD.sum() > 0) {
    Mu = M_.selfadjointView<Eigen::Upper>() * u + u;
    num = utils::selectFromIndicator(Mu, idxD);
    den = utils::selectFromIndicator(Cbu, idxD);
    d = (num.array() / den.array()).minCoeff();    
  
    #ifdef DEBUG_OPTIM
      std::cout << "clipper: u0: ";
      for (int ii=0; ii<u0.size(); ++ii) {std::cout << u0(ii) << " ";} std::cout << std::endl;

      std::cout << "clipper: u: ";
      for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";} std::cout << std::endl;

      std::cout << "clipper: initial d: " << d << std::endl;
    #endif
  }

  #ifdef DEBUG_TIMING
      const auto t2 = std::chrono::high_resolution_clock::now(); // timer
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t_start);
      auto elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: u0, d set up: " << elaps << std::endl;
  #endif

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: params_.tol_u: " << params_.tol_u << std::endl;
    std::cout << "clipper: params_.tol_F: " << params_.tol_F << std::endl;
    std::cout << "clipper: params_.eps: " << params_.eps << std::endl;
    std::cout << "clipper: params_.maxlsiters: " << params_.maxlsiters << std::endl;    
  #endif

  //NOTE: we assume matrices M_ & C_ are identical & binary

  // Orthogonal projected gradient ascent 
  double F = 0; // objective value
  size_t i;
  size_t j;
  size_t k; // iteration counters
  double alpha = 1; // initialize stepsize for line search
  for (i=0; i<params_.maxoliters; ++i) { // outerloop
    
    uold = u; // store last solution
    gradF = (1 + d) * u - d * ones * u.sum() + M_.selfadjointView<Eigen::Upper>() * u * (1 + d);
    F = u.dot(gradF); // current objective value
    gradF = gradF - F * u; // project gradient onto S^n tangent bundle (note here F = u.dot(gradF))

    for (j=0; j<params_.maxiniters; ++j) { // innerloop                        
      
      // Identify an aggressive initial step size based on which elements
      // the gradient indicates can be penalized (unless already negative)      
      const Eigen::VectorXi idxA = ((gradF.array() < -params_.eps) && ( u.array() > params_.eps)).cast<int>();
      if (idxA.sum() > 0) {
        num = utils::selectFromIndicator(u, idxA);
        den = utils::selectFromIndicator(gradF, idxA);
        alpha = (num.array() / den.array()).abs().minCoeff(); // smallest alpha that leads to a step that hits the positive orthant boundaries (for indices where gradF<0 and u>0)  
      } else {
        alpha = std::pow((1.0/params_.beta), 3) / gradF.norm(); // 3 steps of backtracking line search if step is too large before gradu becomes norm 1
      }

      // Backtracking line search on gradient ascent
      double Fnew = 0;
      double deltaF = 0;
      for (k=0; k<params_.maxlsiters; ++k) {
        unew = u + alpha * gradF;                     // gradient step
        unew = unew.cwiseMax(0);                      // project onto positive orthant
        unew.normalize();                             // project onto S^n
        gradFnew = (1 + d) * unew // because M/C is missing identity on diagonal
                    - d * ones * unew.sum()
                    + M_.selfadjointView<Eigen::Upper>() * unew * (1 + d);
        Fnew = unew.dot(gradFnew);                    // new objective value after step
        gradFnew = gradFnew - Fnew * unew;            // project gradient

        deltaF = Fnew - F;                            // change in objective value

        if (deltaF < -params_.eps) {
          // objective value decreased---we need to backtrack, so reduce step size
          alpha = alpha * params_.beta;
        } else {
          break; // obj value increased, stop line search
        }
        

        #ifdef DEBUG_OPTIM
          if (k==params_.maxlsiters-1) {std::cout << "clipper: reached max line search itr." << std::endl;}
        #endif
      } // end line search
      const double deltau = (unew.array()-u.array()).abs().maxCoeff(); // l1 norm
      // update values
      F = Fnew;
      u = unew;
      gradF = gradFnew;

      // check if desired accuracy has been reached by gradient ascent 
      if (deltau < params_.tol_u || std::abs(deltaF) < params_.tol_F) break;      
    } // end innerloop

    #ifdef DEBUG_OPTIM
      std::cout << "clipper: j loop: " << j << std::endl; 
      std::cout << "clipper: gradf: ";
      for (int ii=0; ii<gradF.size(); ++ii) {std::cout << gradF(ii) << " ";} std::cout << std::endl;     
    #endif

    // exit outer loop once changes in u are small
    const double deltau_elem = (u.array() - uold.array()).abs().maxCoeff();

    // Increase d
    Cbu = ones * u.sum() - C_.selfadjointView<Eigen::Upper>() * u - u;
    const Eigen::VectorXi idxD = ((Cbu.array() > params_.eps) && (u.array() > params_.eps)).cast<int>();
    if (idxD.sum() > 0) {
      Mu = M_.selfadjointView<Eigen::Upper>() * u + u;
      num = utils::selectFromIndicator(Mu, idxD);
      den = utils::selectFromIndicator(Cbu, idxD);
      const double deltad = (num.array() / den.array()).abs().minCoeff();

      d += deltad;

      #ifdef DEBUG_OPTIM
        std::cout << "clipper: u: ";
        for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";}
        std::cout << std::endl;
      #endif
    } else {
      break;
    }
    #ifdef DEBUG_OPTIM
      std::cout << "clipper: increased d: " << d << std::endl;
    #endif
  } // end outerloop

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: i loop: " << i << std::endl;
  #endif
  #ifdef DEBUG_TIMING
    const auto t3 = std::chrono::high_resolution_clock::now(); // timer
    dur = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
    elaps = static_cast<double>(dur.count()) / 1e6;
    std::cout << "clipper time: optimization: " << elaps << std::endl;
  #endif
  #ifdef DEBUG_OPTIM
    std::cout << "clipper: u final: ";
    for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";} std::cout << std::endl;
  #endif

  // Generate output
  std::vector<long> nodes; // node indices of rounded u vector

  // pick a rounding threshold between min and max element
  const double rounding_thresh = ( u.array().maxCoeff() - u.array().minCoeff() ) / 2;
  #ifdef DEBUG_OPTIM
    std::cout << "clipper: rounding_thresh: " << rounding_thresh << std::endl;
  #endif

  nodes = utils::findIndicesWhereAboveThreshold(u, rounding_thresh);

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: nodes above eps thresh: ";
    for (long node : nodes) {std::cout << node << " ";} std::cout << std::endl;
  #endif
  #ifdef DEBUG_TIMING
      const auto t4 = std::chrono::high_resolution_clock::now(); // timer
      dur = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
      elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: rounding output: " << elaps << std::endl;
  #endif

  const auto t_end = std::chrono::high_resolution_clock::now();
  const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
  const double elapsed = static_cast<double>(duration.count()) / 1e9;

  // set solution
  soln_.t = elapsed;
  soln_.ifinal = i;
  std::swap(soln_.nodes, nodes);
  soln_.u0 = u0;
  soln_.u.swap(u);
  soln_.score = F;
}

////////////////////////////////////////////////////////////////////////////////////
// find maximal clique for a binray adjacency matrix (i.e., unweighted graph)
// the old version (identical to matlab) 
void CLIPPER::findMaximalClique_old(const Eigen::VectorXd& u0)
{
  const auto t_start = std::chrono::high_resolution_clock::now();
    
  // Initialization
  const size_t n = _M.cols();
  const Eigen::MatrixXd C = _C;

  // Zero out any entry corresponding to an active constraint
  const Eigen::MatrixXd M = _M.cwiseProduct(C);

  // Binary complement of constraint matrix
  const Eigen::MatrixXd Cb = Eigen::MatrixXd::Ones(n,n) - C;

  // one step of power method to have a good scaling of u
  Eigen::VectorXd u = M * u0;
  u.normalize();

  // initial value of d
  double d = 0; // zero if there are no active constraints
  Eigen::MatrixXd Cbu = Cb * u;
  Eigen::MatrixXd Mu = M * u;
  auto num = Eigen::VectorXd(n);
  auto den = Eigen::VectorXd(n);
  Eigen::VectorXi idxD = ((Cbu.array()>params_.eps) && (u.array()>params_.eps)).cast<int>();
  if (idxD.sum() > 0) {    
    num = utils::selectFromIndicator(Mu, idxD);
    den = utils::selectFromIndicator(Cbu, idxD);
    d = (num.array() / den.array()).minCoeff();
  }

  auto Md = Eigen::MatrixXd(M.rows(), M.cols());
  Md = M - d*Cb;

  #ifdef DEBUG_TIMING
      const auto t2 = std::chrono::high_resolution_clock::now(); // timer
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t_start);
      auto elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: u0, d set up: " << elaps << std::endl;
  #endif
  #ifdef DEBUG_OPTIM
    std::cout << "clipper: u0: ";
    for (int ii=0; ii<u0.size(); ++ii) {std::cout << u0(ii) << " ";} std::cout << std::endl;
    std::cout << "clipper: u: ";
    for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";} std::cout << std::endl;
    std::cout << "clipper: initial d: " << d << std::endl;

    std::cout << "clipper: params_.tol_u: " << params_.tol_u << std::endl;
    std::cout << "clipper: params_.tol_F: " << params_.tol_F << std::endl;
    std::cout << "clipper: params_.eps: " << params_.eps << std::endl;
    std::cout << "clipper: params_.maxlsiters: " << params_.maxlsiters << std::endl;    
  #endif
  
  // initialize memory
  auto gradF = Eigen::VectorXd(n);
  auto gradFnew = Eigen::VectorXd(n);
  auto unew = Eigen::VectorXd(n);

  // Orthogonal projected gradient ascent 
  double F = 0; // objective value
  size_t i;
  size_t k; // iteration counters
  size_t jsum=0;
  size_t ksum=0; // total iteration counters
  double alpha = 1; // initialize
  for (i=0; i<params_.maxoliters; ++i) { // outerloop
    gradF = Md * u;
    F = u.dot(gradF); // current objective value
    gradF = gradF - F * u; // orthogonal projection of gradient onto tangent plane to S^n at u

    // Orthogonal projected gradient ascent
    for (size_t j=0; j<params_.maxiniters; ++j) { // innerloop
      // Backtracking line search on gradient ascent
      double Fnew = 0;
      double deltaF = 0;
      for (k=0; k<params_.maxlsiters; ++k) { // line search
        unew = u + alpha * gradF;                     // gradient step
        unew = unew.cwiseMax(0);                      // project onto positive orthant
        unew.normalize();                             // project onto S^n
        gradFnew = Md * unew;                         // new gradient 
        Fnew = unew.dot(gradFnew);                    // new objective value after step
        
        gradFnew = gradFnew - Fnew * unew;            // project gradient onto S^n
        
        deltaF = Fnew - F;                            // change in objective value

        // armijo condition
        bool armijo_cond =  ( deltaF >= params_.sigma * gradF.dot(unew-u) );
        if (!armijo_cond) {        
          alpha = alpha * params_.beta; // backtrack & reduce step size
        } else {
          alpha = alpha / sqrt(params_.beta); // increase step size
          break; // stop line search
        }
        ksum++;

        #ifdef DEBUG_OPTIM
          if (k==params_.maxlsiters-1) {std::cout << "clipper: reached max line search itr." << std::endl;}
        #endif
      } // line search
      
      // chech alpha after line search
      if (alpha < params_.minalpha) { // alpha decreased too much
        alpha = params_.minalpha;
        #ifdef DEBUG_OPTIM
          std::cout << "lower bounded alpha to min threshold"<< std::endl;
        #endif
      } else if (alpha > params_.maxalpha) { // alpha increased too much
        alpha = params_.maxalpha;
        #ifdef DEBUG_OPTIM
          std::cout << "upper bounded alpha to max threshold"<< std::endl;
        #endif        
      }
      
      const double deltau = (unew-u).norm(); // change in vector u

      // update values
      u = unew;
      F = Fnew;      
      gradF = gradFnew;

      // check if desired accuracy has been reached by gradient ascent 
      if (deltau < params_.tol_u || std::abs(deltaF) < params_.tol_F) break;
      jsum++;
    } //innerloop

    // Increase d
    Cbu = Cb * u;        
    idxD = ((Cbu.array()>params_.eps) && (u.array()>params_.eps)).cast<int>();
    if (idxD.sum() > 0) {      
      Mu = M * u; 
      num = utils::selectFromIndicator(Mu, idxD);
      den = utils::selectFromIndicator(Cbu, idxD);
      const double deltad = (num.array() / den.array()).abs().minCoeff();

      d += deltad;
      Md = M - d*Cb; // update Md
    } else {
      break;
    }
  } // outerloop

  #ifdef DEBUG_TIMING
      const auto t3 = std::chrono::high_resolution_clock::now(); // timer
      dur = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
      elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: optimization: " << elaps << std::endl;
      std::cout << "d: " << d << std::endl;
      std::cout << "ksum: " << ksum << std::endl;
      std::cout << "jsum: " << jsum << std::endl;
      std::cout << "i: " << i << std::endl;      
  #endif
  #ifdef DEBUG_OPTIM
    std::cout << "clipper: u final: ";
    for (int ii=0; ii<u.size(); ++ii) {std::cout << u(ii) << " ";} std::cout << std::endl;
  #endif

  // Generate output
  std::vector<long> nodes; // node indices of rounded u vector

  // pick a rounding threshold between min and max element
  const double rounding_thresh = params_.eps;
  #ifdef DEBUG_OPTIM
    std::cout << "clipper: rounding_thresh: " << rounding_thresh << std::endl;
  #endif

  nodes = utils::findIndicesWhereAboveThreshold(u, rounding_thresh);

  #ifdef DEBUG_OPTIM
    std::cout << "clipper: nodes above eps thresh: ";
    for (long i : nodes) {std::cout << i << " ";} std::cout << std::endl;
  #endif
  #ifdef DEBUG_TIMING
      const auto t4 = std::chrono::high_resolution_clock::now(); // timer
      dur = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
      elaps = static_cast<double>(dur.count()) / 1e6;
      std::cout << "clipper time: rounding output: " << elaps << std::endl;
  #endif

  const auto t_end = std::chrono::high_resolution_clock::now();
  const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
  const double elapsed = static_cast<double>(duration.count()) / 1e6;

  // set solution
  soln_.t = elapsed;
  soln_.ifinal = i;
  std::swap(soln_.nodes, nodes);
  soln_.u0 = u0;
  soln_.u.swap(u);
  soln_.score = F;
}

} // ns clipper
