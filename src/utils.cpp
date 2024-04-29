
#include "clipperplus/utils.h"

namespace clipperplus::utils
{

Eigen::VectorXd selectFromIndicator(
    const Eigen::VectorXd &x,
    const Eigen::VectorXi &ind)
{
    Eigen::VectorXd y(ind.sum());
    size_t idx = 0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (ind[i])
        {
            y[idx] = x[i];
            idx++;
        }
    }
    return y;
}

std::vector<long> findIndicesWhereAboveThreshold(
    const Eigen::VectorXd& x,
    double thr
) {
  std::vector<long> indices;
  indices.reserve(x.size());
  for (size_t i=0; i<x.rows(); ++i) {
    if (x(i) > thr) indices.push_back(i);
  }
  return indices;
}


}