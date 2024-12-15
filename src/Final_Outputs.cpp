//
//  Final_Outputs.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/11/24.
//

#include "Final_Outputs.hpp"
#include "IOD.hpp"
#include "Arrival_Orbit.hpp" 
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include "Plot_Earth_Moon_Rotating.hpp"
#include "Recreating_Transfer.hpp"
#include "Orbit_Plot.hpp"

void Create_Outputs(IOD& state1, Arrival_Orbit& arrival_orbit, double LU, double TU, double mu, double x_departure, double y_departure, double z_departure, double xdot_departure, double ydot_departure, double zdot_departure, double T1, double delta_vx, double delta_vy, double delta_vz, double T2_elapsed, bool PrintMetaData){
    using namespace std;
    typedef std::vector< double > state_type;

    std::ofstream outputFile("/Users/jonathonhope/Desktop/Cpp 424/Final_Project_Outputs/Test6", std::ios::app);

    if (!outputFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing." << std::endl;
            return;
        }
    
    if(PrintMetaData){
        string Orbit_type = "Out-of-Plane Resonant";
        outputFile << "Meta_Start" << endl << endl << "Arrival Orbit" << '\t' << '\t' << "= "<< arrival_orbit.Orbit_Type << ", orbit number:" << arrival_orbit.Arrival_orbit_number << endl << "Epsilon_km" << '\t' << '\t' << "= 70" << endl << "LU_km" << '\t' << '\t' << "= " << LU << endl << "TU_s" << '\t' << '\t' << "= " << TU << endl << "mu_LU_TU" << '\t' << '\t' << "= " << std::fixed << std::setprecision(12) << mu << endl << "Departure Position:" << '\t' << '\t' << "["<< state1.r2_vec[0] << '\t' << state1.r2_vec[1] << '\t' << state1.r2_vec[2] <<"]" << endl << "Departure Velocity:" << '\t' << '\t' << "["<< state1.v2_vec[0] << '\t' << state1.v2_vec[1] << '\t' << state1.v2_vec[2] <<"]" <<"Meta_End" << endl << endl;
        
        outputFile << "x_departure_LU" << '\t' << "y_departure_LU" << '\t' << "z_departure_LU" << '\t' << "xdot_departure_LU_TU" << '\t' << "ydot_departure_LU_TU" << '\t' << "zdot_departure_LU_TU" << '\t' << "T1_hours" << '\t' << "delta_v2_x_LU_TU" << '\t' << "delta_v2_y_LU_TU"<< '\t' << "delta_v2_z_LU_TU" << '\t' << "T2_hours" << endl;
        
        
        
    }
        
        
        
        cout << std::fixed << std::setprecision(14) << x_departure << '\t' << y_departure << '\t' <<  z_departure << '\t' << xdot_departure << '\t' << ydot_departure << '\t' << zdot_departure << '\t' << T1 << '\t' << delta_vx << '\t' << delta_vy << '\t' << delta_vz << '\t' << T2_elapsed << '\n';
        outputFile << std::fixed << std::setprecision(14) << x_departure << '\t' << y_departure << '\t' <<  z_departure << '\t' << xdot_departure << '\t' << ydot_departure << '\t' << zdot_departure << '\t' << T1 << '\t' << delta_vx << '\t' << delta_vy << '\t' << delta_vz << '\t' << T2_elapsed << '\n';
    
    
}



