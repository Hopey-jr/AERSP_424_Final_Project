//
//  PseudoInverse.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "PseudoInverse.hpp"
#include <Eigen/Eigenvalues>

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd &mat, double tolerance) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance_val = tolerance * std::max(mat.cols(), mat.rows()) * svd.singularValues().array().abs()(0);
    Eigen::ArrayXd singularValuesInv = svd.singularValues().array().abs().inverse();
    singularValuesInv = (svd.singularValues().array().abs() > tolerance_val).select(singularValuesInv, 0);
    return svd.matrixV() * singularValuesInv.matrix().asDiagonal() * svd.matrixU().adjoint();
}
