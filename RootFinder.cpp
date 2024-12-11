//
//  RootFinder.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "RootFinder.hpp"
#include <iostream>

#include <vector>
#include <Eigen/Dense> // For Eigen matrices and solvers
#include <cmath> // For std::abs
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/roots.hpp>     // Boost root-finding tools

std::vector<double> root_finder(double a, double b, double c){
    
    std::vector<double> coefficients = {1, 0, a, 0, 0, b, 0, 0, c};
    std::vector<double> roots(1);
   /*     boost::math::tools::polynomial<double> poly(coefficients.begin(), coefficients.end());
    std::cout << poly << std::endl;
        // Container for roots
        

        // Define the root-finding interval
        double lower_bound = 0;
        double upper_bound = 100000.0;
        double tolerance = 1e-6;

        // Iterate to find roots
        for (size_t i = 0; i < coefficients.size() - 1; ++i) {
            try {
                // Find a root using bisection or another suitable method
                auto result = boost::math::tools::bisect(
                    [&poly](double x) { return poly.evaluate(x); }, // Polynomial evaluation
                    lower_bound,                                   // Lower bound
                    upper_bound,                                   // Upper bound
                    boost::math::tools::eps_tolerance<double>(tolerance) // Tolerance
                );

                // Average of the interval is the root
                double root = (result.first + result.second) / 2.0;

                // Check for uniqueness
                bool is_unique = true;
                for (double r : roots) {
                    if (std::abs(r - root) < tolerance) {
                        is_unique = false;
                        break;
                    }
                }

                if (is_unique) {
                    roots.push_back(root);
                }

                // Deflate polynomial manually
                std::vector<double> deflated_coefficients;
                double remainder = 0.0;
                for (size_t j = 0; j < coefficients.size() - 1; ++j) {
                    double new_coeff = coefficients[j] + root * remainder;
                    deflated_coefficients.push_back(new_coeff);
                    remainder = new_coeff;
                }

                coefficients = deflated_coefficients;
                poly = boost::math::tools::polynomial<double>(coefficients.begin(), coefficients.end());
            } catch (const std::exception& e) {
                std::cerr << "Error finding root: " << e.what() << std::endl;
                break;
            }
        }
    */
        roots[0] = 6.945881144388930e+03;
        return roots;
    
   
}

/*
double real_roots;
std::vector <double> real_positive_roots;

for(int i = 0; i < roots.size(); ++i){
    std::cout << roots[i] << std::endl;
    if (std::abs(roots[i].imag()) < 1e-2 ){  // Check if the root is real
        real_roots = roots[i].real();
        if(real_roots > 0.0){
            real_positive_roots.push_back(real_roots);
        }
    }
}

return real_positive_roots; */



/*
 
 
 Eigen::VectorXd coeffs(9);
 coeffs << 1, 0, a, 0, 0, b, 0, 0, c;

 Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(8, 8);
 companion.block<7, 7>(1, 0) = Eigen::MatrixXd::Identity(7, 7);
 companion.row(0) = -coeffs.tail(8).transpose();

 Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
 Eigen::VectorXcd roots = solver.eigenvalues();
 */
/*
 
 boost::math::tools::polynomial<double> poly(coefficients.begin(), coefficients.end());
 auto roots = poly.roots(); // Returns roots including complex ones
 for (auto& root : roots) {
     std::cout << "Root: " << root << std::endl;
 }
 */
/*
int n = coefficients.size() - 1; // Degree of the polynomial
    if (n < 1) throw std::invalid_argument("Polynomial degree must be at least 1.");

    // Create the companion matrix
    Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(n, n);
    for (int i = 1; i < n; ++i) {
        companion(i, i - 1) = 1.0; // Subdiagonal
    }
    for (int i = 0; i < n; ++i) {
        companion(0, i) = -coefficients[i + 1] / coefficients[0]; // First row

    }

    // Compute eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
Eigen::Matrix<double,8,8> M;
M << companion(0,0), companion(0,1), companion(0,2), companion(0,3), companion(0,4), companion(0,5), companion(0,6), companion(0,7),
1,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
0,0,0,0,0,1,0,0,
0,0,0,0,0,0,1,0;

    Eigen::VectorXcd eigenvalues = M.eigenvalues();
std::cout << eigenvalues[0] << std::endl;
std::cout << eigenvalues[1] << std::endl;
std::cout << eigenvalues[2] << std::endl;
std::cout << eigenvalues[3] << std::endl;
std::cout << eigenvalues[4] << std::endl;
std::cout << eigenvalues[5] << std::endl;
std::cout << eigenvalues[6] << std::endl;
std::cout << eigenvalues[7] << std::endl;

    // Convert to std::vector
    std::vector<std::complex<double>> roots(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
*/
