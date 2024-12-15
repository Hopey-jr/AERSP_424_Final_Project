//  Description:
//               The purpose of this code is to propagate the equations of
//               motion from t0 to tf and is supposed to be used with a numerical integrator.


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

//function call:
// Refer to the specific ode solver for the function call

#include "EOM_function.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>


typedef std::vector< double > state_type;
namespace ublas = boost::numeric::ublas;

void state_prop( const state_type &w , state_type &dwdt , const double /* t */ )
{
     double mu = 1.215058560962404E-2;

    dwdt[0] = w[3];
    dwdt[1] = w[4];
    dwdt[2] = w[5];
 
    double P1;
    double P2;
    double Uy;
    double Ux;
     P1=pow(pow(w[0]+mu,2)+pow(w[1],2)+pow(w[2],2),.5);
     P2=pow(pow(w[0]-1+mu,2)+pow(w[1],2)+pow(w[2],2),.5);
     Ux=w[0]-((1-mu)*(w[0]+mu))/pow(P1,3)-(mu*(w[0]-1+mu))/pow(P2,3);
     Uy=w[1]-((1-mu)*w[1])/pow(P1,3)-mu*w[1]/pow(P2,3);
    dwdt[3] =  2*w[4]+Ux;
    dwdt[4] =  -2*w[3]+Uy;
    dwdt[5] = -(1-mu)*w[2]/pow(P1,3)-mu*w[2]/pow(P2,3);
    
    
    
}
