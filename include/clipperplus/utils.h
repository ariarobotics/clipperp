#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>

namespace clipperplus::utils
{

/**
 * @brief      Select the elements of a vector x given an indicator vector.
 *
 * @param[in]  x     Vector to select elements of
 * @param[in]  ind   The indicator vector
 *
 * @return     Vector of selected elements, with size <= x.size
 */
Eigen::VectorXd selectFromIndicator(
    const Eigen::VectorXd& x,
    const Eigen::VectorXi& ind);



std::vector<long> findIndicesWhereAboveThreshold(
    const Eigen::VectorXd& x,
    double thr
);
}
