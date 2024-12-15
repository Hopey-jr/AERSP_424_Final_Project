//Description:
//            The purpose of this function is to find the pseudo-inverse of a non-square matrix using the Eigen library

//Input:
//       mat - The matrix that the user wants to invert
//       tolerance - The tolerance determined to be valid to produce a valid inverse

//Output:
//       pinv - The inverse of the non-square matrix

//Local:



//Function Call:
//              Eigen::MatrixXd pinv = pseudoInverse(Augmented_STM);



#include "PseudoInverse.hpp"
#include <Eigen/Eigenvalues>

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd &mat, double tolerance) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance_val = tolerance * std::max(mat.cols(), mat.rows()) * svd.singularValues().array().abs()(0);
    Eigen::ArrayXd singularValuesInv = svd.singularValues().array().abs().inverse();
    singularValuesInv = (svd.singularValues().array().abs() > tolerance_val).select(singularValuesInv, 0);
    return svd.matrixV() * singularValuesInv.matrix().asDiagonal() * svd.matrixU().adjoint();
}
