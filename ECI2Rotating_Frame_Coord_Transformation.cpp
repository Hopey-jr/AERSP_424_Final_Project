//
//  ECI2Rotating_Frame_Coord_Transformation.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#include "ECI2Rotating_Frame_Coord_Transformation.hpp"

#include <iostream>
#include <vector>
#include <cmath>

std::vector<double>ECI2Rotating_Frame_Rotation(double x_ECI, double y_ECI ,double z_ECI, double xdot_ECI, double ydot_ECI ,double zdot_ECI, double t){
    
   double LU_JPL = 389703;
   double TU_JPL = 382981;
   double mu_E = pow(LU_JPL,3)/pow(TU_JPL,2);
   double LU = 384400;
   double TU = pow((pow(LU,3)/mu_E),0.5);
   double mu = 1.215058560962404E-2;
    
    //[ non_dimesionalizing the inputs
    double tau = t/TU;
     x_ECI = x_ECI/LU;
     y_ECI = y_ECI/LU;
     z_ECI = z_ECI/LU;
       
     xdot_ECI = xdot_ECI/LU*TU;
     ydot_ECI = ydot_ECI/LU*TU;
     zdot_ECI = zdot_ECI/LU*TU;
       
       //]
    
    //[ The coordinates are rotated into the rotating frame and then translated to the Earth at (-mu, 0, 0)
    double x_rot = x_ECI*cos(tau)+y_ECI*sin(tau) - mu;
    double y_rot = -x_ECI*sin(tau)+y_ECI*cos(tau);
    double z_rot = z_ECI;
       
    double xdot_rot = xdot_ECI*cos(tau)-x_ECI*sin(tau)+ydot_ECI*sin(tau)+y_ECI*cos(tau);
    double ydot_rot = -xdot_ECI*sin(tau)-x_ECI*cos(tau)+ydot_ECI*cos(tau)-y_ECI*sin(tau);
    double zdot_rot = zdot_ECI;
       //]
    
    std::vector<double> r_rot(6);
    r_rot[0] = x_rot;
    r_rot[1] = y_rot;
    r_rot[2] = z_rot;
    r_rot[3] = xdot_rot;
    r_rot[4] = ydot_rot;
    r_rot[5] = zdot_rot;
    return r_rot;
}
