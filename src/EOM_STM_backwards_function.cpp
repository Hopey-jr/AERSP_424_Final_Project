//  Description:
//               The purpose of this code is to propagate the equations of
//               motion and the state transition matrix from t0 to Tf and is
//               supposed to be used with a numerical integrator. This code, as
//               opposed to 'EOM_STM_function' is used when the state
//               transition matrix needs to be propagated in backwards time.


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

#include "EOM_STM_backwards_function.hpp"
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

void state_STM_prop_backwards(const state_type &w, state_type &dwdt, const double /* t */)
{
    double mu = 1.215058560962404E-2;

    dwdt[0] = -w[3];
    dwdt[1] = -w[4];
    dwdt[2] = -w[5];

    double P1;
    double P2;
    double Uy;
    double Ux;
    P1 = pow(pow(w[0] + mu, 2) + pow(w[1], 2) + pow(w[2], 2), .5);
    P2 = pow(pow(w[0] - 1 + mu, 2) + pow(w[1], 2) + pow(w[2], 2), .5);
    Ux = w[0] - ((1 - mu) * (w[0] + mu)) / pow(P1, 3) - (mu * (w[0] - 1 + mu)) / pow(P2, 3);
    Uy = w[1] - ((1 - mu) * w[1]) / pow(P1, 3) - mu * w[1] / pow(P2, 3);
    dwdt[3] = -1 * (2 * w[4] + Ux);
    dwdt[4] = -1 * (-2 * w[3] + Uy);
    dwdt[5] = -1 * (-(1 - mu) * w[2] / pow(P1, 3) - mu * w[2] / pow(P2, 3));

    double Uxx = 1 - (1 - mu) / pow(P1, 3) - mu / pow(P2, 3) + 3 * (1 - mu) * pow(w[0] + mu, 2) / pow(P1, 5) + 3 * mu * pow(w[0] - 1 + mu, 2) / pow(P2, 5);
    double Uyy = 1 - (1 - mu) / pow(P1, 3) - mu / pow(P2, 3) + 3 * (1 - mu) * pow(w[1], 2) / pow(P1, 5) + 3 * mu * pow(w[1], 2) / pow(P2, 5);
    double Uzz = -(1 - mu) / pow(P1, 3) - mu / pow(P2, 3) + 3 * (1 - mu) * pow(w[2], 2) / pow(P1, 5) + 3 * mu * pow(w[2], 2) / pow(P2, 5);

    double Uxy = 3 * (1 - mu) * (w[0] + mu) * w[1] / pow(P1, 5) + 3 * mu * (w[0] - 1 + mu) * w[1] / pow(P2, 5);
    double Uxz = 3 * (1 - mu) * (w[0] + mu) * w[2] / pow(P1, 5) + 3 * mu * (w[0] - 1 + mu) * w[2] / pow(P2, 5);
    double Uyz = 3 * (1 - mu) * w[1] * w[2] / pow(P1, 5) + 3 * mu * w[1] * w[2] / pow(P2, 5);

    ublas::matrix<double> STM(6, 6);
    ublas::matrix<double> STM_dot(6, 6);
    ublas::matrix<double> A(6, 6);

    A(0, 0) = 0;
    A(0, 1) = 0;
    A(0, 2) = 0;
    A(0, 3) = 1;
    A(0, 4) = 0;
    A(0, 5) = 0;

    A(1, 0) = 0;
    A(1, 1) = 0;
    A(1, 2) = 0;
    A(1, 3) = 0;
    A(1, 4) = 1;
    A(1, 5) = 0;

    A(2, 0) = 0;
    A(2, 1) = 0;
    A(2, 2) = 0;
    A(2, 3) = 0;
    A(2, 4) = 0;
    A(2, 5) = 1;

    A(3, 0) = Uxx;
    A(3, 1) = Uxy;
    A(3, 2) = Uxz;
    A(3, 3) = 0;
    A(3, 4) = 2;
    A(3, 5) = 0;

    A(4, 0) = Uxy;
    A(4, 1) = Uyy;
    A(4, 2) = Uyz;
    A(4, 3) = -2;
    A(4, 4) = 0;
    A(4, 5) = 0;

    A(5, 0) = Uxz;
    A(5, 1) = Uyz;
    A(5, 2) = Uzz;
    A(5, 3) = 0;
    A(5, 4) = 0;
    A(5, 5) = 0;

    STM(0, 0) = w[6];
    STM(0, 1) = w[7];
    STM(0, 2) = w[8];
    STM(0, 3) = w[9];
    STM(0, 4) = w[10];
    STM(0, 5) = w[11];

    STM(1, 0) = w[12];
    STM(1, 1) = w[13];
    STM(1, 2) = w[14];
    STM(1, 3) = w[15];
    STM(1, 4) = w[16];
    STM(1, 5) = w[17];

    STM(2, 0) = w[18];
    STM(2, 1) = w[19];
    STM(2, 2) = w[20];
    STM(2, 3) = w[21];
    STM(2, 4) = w[22];
    STM(2, 5) = w[23];

    STM(3, 0) = w[24];
    STM(3, 1) = w[25];
    STM(3, 2) = w[26];
    STM(3, 3) = w[27];
    STM(3, 4) = w[28];
    STM(3, 5) = w[29];

    STM(4, 0) = w[30];
    STM(4, 1) = w[31];
    STM(4, 2) = w[32];
    STM(4, 3) = w[33];
    STM(4, 4) = w[34];
    STM(4, 5) = w[35];

    STM(5, 0) = w[36];
    STM(5, 1) = w[37];
    STM(5, 2) = w[38];
    STM(5, 3) = w[39];
    STM(5, 4) = w[40];
    STM(5, 5) = w[41];

    STM_dot = ublas::prod(A, STM);

    dwdt[6] = -STM_dot(0, 0);
    dwdt[7] = -STM_dot(0, 1);
    dwdt[8] = -STM_dot(0, 2);
    dwdt[9] = -STM_dot(0, 3);
    dwdt[10] = -STM_dot(0, 4);
    dwdt[11] = -STM_dot(0, 5);

    dwdt[12] = -STM_dot(1, 0);
    dwdt[13] = -STM_dot(1, 1);
    dwdt[14] = -STM_dot(1, 2);
    dwdt[15] = -STM_dot(1, 3);
    dwdt[16] = -STM_dot(1, 4);
    dwdt[17] = -STM_dot(1, 5);

    dwdt[18] = -STM_dot(2, 0);
    dwdt[19] = -STM_dot(2, 1);
    dwdt[20] = -STM_dot(2, 2);
    dwdt[21] = -STM_dot(2, 3);
    dwdt[22] = -STM_dot(2, 4);
    dwdt[23] = -STM_dot(2, 5);

    dwdt[24] = -STM_dot(3, 0);
    dwdt[25] = -STM_dot(3, 1);
    dwdt[26] = -STM_dot(3, 2);
    dwdt[27] = -STM_dot(3, 3);
    dwdt[28] = -STM_dot(3, 4);
    dwdt[29] = -STM_dot(3, 5);

    dwdt[30] = -STM_dot(4, 0);
    dwdt[31] = -STM_dot(4, 1);
    dwdt[32] = -STM_dot(4, 2);
    dwdt[33] = -STM_dot(4, 3);
    dwdt[34] = -STM_dot(4, 4);
    dwdt[35] = -STM_dot(4, 5);

    dwdt[36] = -STM_dot(5, 0);
    dwdt[37] = -STM_dot(5, 1);
    dwdt[38] = -STM_dot(5, 2);
    dwdt[39] = -STM_dot(5, 3);
    dwdt[40] = -STM_dot(5, 4);
    dwdt[41] = -STM_dot(5, 5);
}
