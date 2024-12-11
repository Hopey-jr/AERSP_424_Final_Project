//
//  EOM_backwards_function.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "EOM_backwards_function.hpp"
#include <vector>


#include <vector>  // For std::vector
#include <cmath>  // For std::pow
typedef std::vector<double> state_type;

void state_prop_backwards(const state_type &w, state_type &dwdt, const double /* t */) {
    double mu = 1.215058560962404E-2;  // Standard gravitational parameter

    // State equations for position derivatives
    dwdt[0] = -w[3];
    dwdt[1] = -w[4];
    dwdt[2] = -w[5];
 
    double P1;
    double P2;
    double Uy;
    double Ux;
    
    // Calculate potential terms
    P1 = std::pow(std::pow(w[0] + mu, 2) + std::pow(w[1], 2) + std::pow(w[2], 2), 0.5);
    P2 = std::pow(std::pow(w[0] - 1 + mu, 2) + std::pow(w[1], 2) + std::pow(w[2], 2), 0.5);
    
    Ux = w[0] - ((1 - mu) * (w[0] + mu)) / std::pow(P1, 3) - (mu * (w[0] - 1 + mu)) / std::pow(P2, 3);
    Uy = w[1] - ((1 - mu) * w[1]) / std::pow(P1, 3) - mu * w[1] / std::pow(P2, 3);
    
    // Derivatives for velocity components
    dwdt[3] = -1 * (2 * w[4] + Ux);
    dwdt[4] = -1 * (-2 * w[3] + Uy);
    dwdt[5] = -1 * (-(1 - mu) * w[2] / std::pow(P1, 3) - mu * w[2] / std::pow(P2, 3));
}


