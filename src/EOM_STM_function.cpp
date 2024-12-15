//  Description:
//               The purpose of this code is to propagate the equations of
//               motion and the state transition matrix from t0 to Tf and is
//               supposed to be used with a numerical integrator. 


//Inputs:
//NOTE: Since this code is supposed to be used in conjunction with an ode
//       solver,  please refer to the specific ode solver to see what inputs
//       are needed. Below is what is generally needed
//       t0 - The initial time of the propagation
//       tf - The final time of the propagation
//       dt - The timestep of the integration between t0 and tf
//       x0 - The initial conditions for the periodic orbit coming from the
//       orbit catalog

//Output:
//       w_vec - The state (position and velocity) and state transition matrix
//       at timesteps from t0 to tf
//       times - The timestamp at each timestep

//Locals:
//       mu - The nondimensional reduced mass of the Earth-Moon system [LU^3/TU^2]
//       w - The state at each timestep
//       dwdt - The time derivative of the state at each timestep
//       p1 - The distance from the Earth to the barrycenter
//       p2 - The distance from the Moon to the barrycenter
//       Ux - The derivative with respect to x of the pseudo potential
//       Uy - The derivative with respect to y of the pseudo potential
//       Uz - The derivative with respect to z of the pseudo potential
//       Uxx - The derivative of Ux with respect to x
//       Uyy - The derivative of Uy with respect to y
//       Uzz - The derivative of Uz with respect to z
//       Uxy - The derivative of Ux with respect to y
//       Uxz - The derivative of Ux with respect to z
//       Uyz - The derivative of Uy with respect to z
//Function call:
// Refer to the specific ode solver for the function call

#include "EOM_STM_function.h"
#include <cmath>

// Implementation of the state_STM_prop function
void state_STM_prop(const state_type &w, state_type &dwdt, const double /* t */) {
    const double mu = 1.215058560962404E-2;

    // Initialize derivatives
    dwdt[0] = w[3];
    dwdt[1] = w[4];
    dwdt[2] = w[5];

    // Compute distances
    const double P1 = std::pow(std::pow(w[0] + mu, 2) + std::pow(w[1], 2) + std::pow(w[2], 2), 0.5);
    const double P2 = std::pow(std::pow(w[0] - 1 + mu, 2) + std::pow(w[1], 2) + std::pow(w[2], 2), 0.5);

    // Compute gravitational accelerations
    const double Ux = w[0] - ((1 - mu) * (w[0] + mu)) / std::pow(P1, 3) - (mu * (w[0] - 1 + mu)) / std::pow(P2, 3);
    const double Uy = w[1] - ((1 - mu) * w[1]) / std::pow(P1, 3) - (mu * w[1]) / std::pow(P2, 3);

    dwdt[3] = 2 * w[4] + Ux;
    dwdt[4] = -2 * w[3] + Uy;
    dwdt[5] = -(1 - mu) * w[2] / std::pow(P1, 3) - mu * w[2] / std::pow(P2, 3);

    // Compute second derivatives of the potential function
    const double Uxx = 1 - (1 - mu) / std::pow(P1, 3) - mu / std::pow(P2, 3)
                       + 3 * (1 - mu) * std::pow(w[0] + mu, 2) / std::pow(P1, 5)
                       + 3 * mu * std::pow(w[0] - 1 + mu, 2) / std::pow(P2, 5);
    const double Uyy = 1 - (1 - mu) / std::pow(P1, 3) - mu / std::pow(P2, 3)
                       + 3 * (1 - mu) * std::pow(w[1], 2) / std::pow(P1, 5)
                       + 3 * mu * std::pow(w[1], 2) / std::pow(P2, 5);
    const double Uzz = -(1 - mu) / std::pow(P1, 3) - mu / std::pow(P2, 3)
                       + 3 * (1 - mu) * std::pow(w[2], 2) / std::pow(P1, 5)
                       + 3 * mu * std::pow(w[2], 2) / std::pow(P2, 5);
    const double Uxy = 3 * (1 - mu) * (w[0] + mu) * w[1] / std::pow(P1, 5)
                       + 3 * mu * (w[0] - 1 + mu) * w[1] / std::pow(P2, 5);
    const double Uxz = 3 * (1 - mu) * (w[0] + mu) * w[2] / std::pow(P1, 5)
                       + 3 * mu * (w[0] - 1 + mu) * w[2] / std::pow(P2, 5);
    const double Uyz = 3 * (1 - mu) * w[1] * w[2] / std::pow(P1, 5)
                       + 3 * mu * w[1] * w[2] / std::pow(P2, 5);

    // Initialize STM and its derivative matrices
    ublas::matrix<double> STM(6, 6);
    ublas::matrix<double> STM_dot(6, 6);
    ublas::matrix<double> A(6, 6);

    // Fill the A matrix
    A(0, 0) = 0;  A(0, 1) = 0;  A(0, 2) = 0;  A(0, 3) = 1;  A(0, 4) = 0;  A(0, 5) = 0;
    A(1, 0) = 0;  A(1, 1) = 0;  A(1, 2) = 0;  A(1, 3) = 0;  A(1, 4) = 1;  A(1, 5) = 0;
    A(2, 0) = 0;  A(2, 1) = 0;  A(2, 2) = 0;  A(2, 3) = 0;  A(2, 4) = 0;  A(2, 5) = 1;
    A(3, 0) = Uxx;  A(3, 1) = Uxy;  A(3, 2) = Uxz;  A(3, 3) = 0;  A(3, 4) = 2;  A(3, 5) = 0;
    A(4, 0) = Uxy;  A(4, 1) = Uyy;  A(4, 2) = Uyz;  A(4, 3) = -2; A(4, 4) = 0;  A(4, 5) = 0;
    A(5, 0) = Uxz;  A(5, 1) = Uyz;  A(5, 2) = Uzz;  A(5, 3) = 0;  A(5, 4) = 0;  A(5, 5) = 0;

    // Fill the STM from the state vector
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            STM(i, j) = w[6 + i * 6 + j];
        }
    }

    // Compute STM_dot
    STM_dot = ublas::prod(A, STM);

    // Fill the derivative vector with STM_dot
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            dwdt[6 + i * 6 + j] = STM_dot(i, j);
        }
    }
}
