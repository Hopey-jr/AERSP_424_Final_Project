//
//  Final_Outputs.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/11/24.
//

#ifndef Final_Outputs_hpp
#define Final_Outputs_hpp

#include "IOD.hpp"          // Include the IOD class header
#include "Arrival_Orbit.hpp" // Include the Arrival_Orbit class header
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>


using namespace std;


void Create_Outputs(IOD& state1, Arrival_Orbit& arrival_orbit, double LU, double TU, double mu, double x_departure, double y_departure, double z_departure, double xdot_departure, double ydot_departure, double zdot_departure, double T1, double delta_vx, double delta_vy, double delta_vz, double T2_elapsed, bool PrintMetaData);
#endif /* Final_Outputs_hpp */
