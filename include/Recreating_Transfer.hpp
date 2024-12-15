//
//  Recreating_Transfer.hpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#ifndef Recreating_Transfer_hpp
#define Recreating_Transfer_hpp

#include "IOD.hpp"          // Include the IOD class header
#include "Arrival_Orbit.hpp" // Include the Arrival_Orbit class header
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>


using namespace std;

typedef std::vector< double > state_type;

std::vector< state_type > Recreate_Transfer(IOD& state1, Arrival_Orbit& arrival_orbit, double LU, double TU, double mu, double x_departure, double y_departure, double z_departure, double xdot_departure, double ydot_departure, double zdot_departure, double T1, double delta_vx, double delta_vy, double delta_vz, double T2_elapsed);
#endif /* Recreating_Transfer_hpp */
