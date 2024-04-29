/**
compute core number of graph vertices, and quickly find a miximal clique.
Code based on Parallel Maximum Clique Algorithm by Ryan A. Rossi (http://ryanrossi.com/pmc)

Author: kaveh fathian (kavehfathian@gmail.com)
 */

#include "clipperplus/clique_optimization.h"

namespace clipperplus
{

    std::vector<long> clique_optimization(
        const Eigen::MatrixXd &_M,
        const Eigen::VectorXd &u0,
        const Params &params)
    {
        const auto t_start = std::chrono::high_resolution_clock::now();

        // Initialization
        const size_t n = _M.cols();
        const Eigen::MatrixXd C = _M;

        // Zero out any entry corresponding to an active constraint
        const Eigen::MatrixXd M = _M.cwiseProduct(C);

        // Binary complement of constraint matrix
        const Eigen::MatrixXd Cb = Eigen::MatrixXd::Ones(n, n) - C;

        // one step of power method to have a good scaling of u
        Eigen::VectorXd u = M * u0;
        u.normalize();

        // initial value of d
        double d = 0; // zero if there are no active constraints
        Eigen::MatrixXd Cbu = Cb * u;
        Eigen::MatrixXd Mu = M * u;
        auto num = Eigen::VectorXd(n);
        auto den = Eigen::VectorXd(n);
        Eigen::VectorXi idxD = ((Cbu.array() > params.eps) && (u.array() > params.eps)).cast<int>();
        if (idxD.sum() > 0)
        {
            num = utils::selectFromIndicator(Mu, idxD);
            den = utils::selectFromIndicator(Cbu, idxD);
            d = (num.array() / den.array()).minCoeff();
        }

        auto Md = Eigen::MatrixXd(M.rows(), M.cols());
        Md = M - d * Cb;

#ifdef DEBUG_TIMING
        const auto t2 = std::chrono::high_resolution_clock::now(); // timer
        auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t_start);
        auto elaps = static_cast<double>(dur.count()) / 1e6;
        std::cout << "clipper time: u0, d set up: " << elaps << std::endl;
#endif
#ifdef DEBUG_OPTIM
        std::cout << "clipper: u0: ";
        for (int ii = 0; ii < u0.size(); ++ii)
        {
            std::cout << u0(ii) << " ";
        }
        std::cout << std::endl;
        std::cout << "clipper: u: ";
        for (int ii = 0; ii < u.size(); ++ii)
        {
            std::cout << u(ii) << " ";
        }
        std::cout << std::endl;
        std::cout << "clipper: initial d: " << d << std::endl;

        std::cout << "clipper: params.tol_u: " << params.tol_u << std::endl;
        std::cout << "clipper: params.tol_F: " << params.tol_F << std::endl;
        std::cout << "clipper: params.eps: " << params.eps << std::endl;
        std::cout << "clipper: params.maxlsiters: " << params.maxlsiters << std::endl;
#endif

        // initialize memory
        auto gradF = Eigen::VectorXd(n);
        auto gradFnew = Eigen::VectorXd(n);
        auto unew = Eigen::VectorXd(n);

        // Orthogonal projected gradient ascent
        double F = 0; // objective value
        size_t i;
        size_t k; // iteration counters
        size_t jsum = 0;
        size_t ksum = 0;  // total iteration counters
        double alpha = 1; // initialize
        for (i = 0; i < params.maxoliters; ++i)
        { // outerloop
            gradF = Md * u;
            F = u.dot(gradF);      // current objective value
            gradF = gradF - F * u; // orthogonal projection of gradient onto tangent plane to S^n at u

            // Orthogonal projected gradient ascent
            for (size_t j = 0; j < params.maxiniters; ++j)
            { // innerloop
                // Backtracking line search on gradient ascent
                double Fnew = 0;
                double deltaF = 0;
                for (k = 0; k < params.maxlsiters; ++k)
                {                              // line search
                    unew = u + alpha * gradF;  // gradient step
                    unew = unew.cwiseMax(0);   // project onto positive orthant
                    unew.normalize();          // project onto S^n
                    gradFnew = Md * unew;      // new gradient
                    Fnew = unew.dot(gradFnew); // new objective value after step

                    gradFnew = gradFnew - Fnew * unew; // project gradient onto S^n

                    deltaF = Fnew - F; // change in objective value

                    // armijo condition
                    bool armijo_cond = (deltaF >= params.sigma * gradF.dot(unew - u));
                    if (!armijo_cond)
                    {
                        alpha = alpha * params.beta; // backtrack & reduce step size
                    }
                    else
                    {
                        alpha = alpha / sqrt(params.beta); // increase step size
                        break;                              // stop line search
                    }
                    ksum++;

#ifdef DEBUG_OPTIM
                    if (k == params.maxlsiters - 1)
                    {
                        std::cout << "clipper: reached max line search itr." << std::endl;
                    }
#endif
                } // line search

                // chech alpha after line search
                if (alpha < params.minalpha)
                { // alpha decreased too much
                    alpha = params.minalpha;
#ifdef DEBUG_OPTIM
                    std::cout << "lower bounded alpha to min threshold" << std::endl;
#endif
                }
                else if (alpha > params.maxalpha)
                { // alpha increased too much
                    alpha = params.maxalpha;
#ifdef DEBUG_OPTIM
                    std::cout << "upper bounded alpha to max threshold" << std::endl;
#endif
                }

                const double deltau = (unew - u).norm(); // change in vector u

                // update values
                u = unew;
                F = Fnew;
                gradF = gradFnew;

                // check if desired accuracy has been reached by gradient ascent
                if (deltau < params.tol_u || std::abs(deltaF) < params.tol_F)
                    break;
                jsum++;
            } // innerloop

            // Increase d
            Cbu = Cb * u;
            idxD = ((Cbu.array() > params.eps) && (u.array() > params.eps)).cast<int>();
            if (idxD.sum() > 0)
            {
                Mu = M * u;
                num = utils::selectFromIndicator(Mu, idxD);
                den = utils::selectFromIndicator(Cbu, idxD);
                const double deltad = (num.array() / den.array()).abs().minCoeff();

                d += deltad;
                Md = M - d * Cb; // update Md
            }
            else
            {
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
        for (int ii = 0; ii < u.size(); ++ii)
        {
            std::cout << u(ii) << " ";
        }
        std::cout << std::endl;
#endif

        // Generate output
        std::vector<long> nodes; // node indices of rounded u vector

        // pick a rounding threshold between min and max element
        const double rounding_thresh = params.eps;
#ifdef DEBUG_OPTIM
        std::cout << "clipper: rounding_thresh: " << rounding_thresh << std::endl;
#endif

        nodes = utils::findIndicesWhereAboveThreshold(u, rounding_thresh);

#ifdef DEBUG_OPTIM
        std::cout << "clipper: nodes above eps thresh: ";
        for (long i : nodes)
        {
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
        const double elapsed = static_cast<double>(duration.count()) / 1e6;

        return nodes;
    }

}
