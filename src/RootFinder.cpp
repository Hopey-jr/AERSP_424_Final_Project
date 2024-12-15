//Description:
//            The purpose of this function is to find the real, positive root of an eighth order polynomial using Newton Raphson iteration
//Input:
//       a - The first coefficient coming from the Gauss function
//       b - The second coefficient coming from the Gauss function
//       c - The third coefficient coming from the Gauss function

//Output:
//       r_new - The distance found as a result of the Newton Raphson algorithm

//Local:
//       r_old - The previous iteration of the distance used in the Newton Raphson algorithm
//       tol - The tolerance needed to hit before the iteration is deemed to converge
//       delta - The difference in the position at each iteration to determine the error

//Function Call:
//       double r2 = root_finder( a,  b,  c);
#include "RootFinder.hpp"
#include <iostream>

#include <vector>
#include <Eigen/Dense> // For Eigen matrices and solvers
#include <cmath> // For std::abs
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/roots.hpp>     // Boost root-finding tools
#include <limits>

using namespace boost::math::tools;

double root_finder(double a, double b, double c){

        const double tol = 1e-8;
       
    double r_old = 6400;
    double r_new = 0;
    double delta = 1;
    //[ Newton Raphson
    
   
    while(delta > tol){
        r_new = r_old - (pow(r_old,8.0)+a*pow(r_old,6.0)+b*pow(r_old,3.0)+c)/(8* pow(r_old,7.0)+6*a* pow(r_old,5.0)+3*b* pow(r_old,2.0));
        delta = abs(r_new-r_old);
        r_old = r_new;
    }
  
        return r_new;
    }
