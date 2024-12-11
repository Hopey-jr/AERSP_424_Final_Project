//
//  Manifold_Negative.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef Manifold_Negative_hpp
#define Manifold_Negative_hpp

#include <vector>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include "PushBackStateAndTime.hpp"

// Define the state type
typedef std::vector<double> state_type;

// Function declaration
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Generate_Negative_Manifolds(const std::vector<double>& Initial_Conditions, double Tp, double mu, int n);
#endif /* Manifold_Negative_hpp */
