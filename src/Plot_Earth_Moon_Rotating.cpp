//
//  Plot_Earth_Moon_Rotating.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#include "Plot_Earth_Moon_Rotating.hpp"

#include <matplot/matplot.h>
#include <cmath>
#include <vector>
#include "ECI2Rotating_Frame_Coord_Transformation.hpp"

void Plot_earth_and_moon(){
    using namespace matplot;

    double pi = M_PI;
   double r_Earth = 6371.0084;
   double r_moon = 1737.4;
   double xdot_ECI = 0.0;
   double ydot_ECI = 0.0;
   double zdot_ECI = 0.0;
    
    const int num_points = 21;
        std::vector<std::vector<double>> X(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Y(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Z(num_points, std::vector<double>(num_points));

        // Define theta and phi ranges
        std::vector<double> theta = matplot::linspace(0, 2 * pi, num_points);
        std::vector<double> phi = matplot::linspace(0, pi, num_points);

        for (size_t i = 0; i < num_points; ++i) {
            for (size_t j = 0; j < num_points; ++j) {
                X[i][j] = cos(theta[i]) * sin(phi[j]);
                Y[i][j] = sin(theta[i]) * sin(phi[j]);
                Z[i][j] = cos(phi[j]);
            }
        }
    
        std::vector<std::vector<double>> X_Earth(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Y_Earth(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Z_Earth(num_points, std::vector<double>(num_points));

        std::vector<std::vector<double>> X_Moon(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Y_Moon(num_points, std::vector<double>(num_points));
        std::vector<std::vector<double>> Z_Moon(num_points, std::vector<double>(num_points));


    
   double time = 0.0;
    for (int j = 0; j < 21; j++){
        for (int i = 0; i < 21; i++){
            std::vector<double> r_rot_Earth = ECI2Rotating_Frame_Rotation(X[i][j]*r_Earth, Y[i][j]*r_Earth, Z[i][j]*r_Earth, xdot_ECI, ydot_ECI, zdot_ECI, time);
    
    
            X_Earth[i][j] = r_rot_Earth[0];
            Y_Earth[i][j] = r_rot_Earth[1];
            Z_Earth[i][j] = r_rot_Earth[2];

    
            std::vector<double> r_rot_Moon = ECI2Rotating_Frame_Rotation(X[i][j]*r_moon, Y[i][j]*r_moon, Z[i][j]*r_moon, xdot_ECI, ydot_ECI, zdot_ECI, time);
            X_Moon[i][j] = r_rot_Moon[0] + 1.0;
            Y_Moon[i][j] = r_rot_Moon[1];
            Z_Moon[i][j] = r_rot_Moon[2];
    
        }
}
    
    surf(X_Earth,Y_Earth,Z_Earth);
    //hold ("on");
    surf(X_Moon,Y_Moon,Z_Moon);
   // hold ("on");
    
    
}
