//
//  PseudoInverse.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef PseudoInverse_hpp
#define PseudoInverse_hpp

#include <Eigen/Dense>

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd &mat, double tolerance = 1e-12);

#endif /* PseudoInverse_hpp */
