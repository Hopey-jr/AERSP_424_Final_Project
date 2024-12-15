//
//  main.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/12/24.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <thread>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <complex>
#include "Arrival_Orbit.hpp"
#include "IOD_Measurements.hpp"
#include "IOD.hpp"
#include "Transfers.hpp"
#include <matplot/matplot.h>
#include "helper.h"
#include "logger_V2.hpp"
#include "Measurement_Input.hpp"
#include "Orbit_Type_Input.hpp"
#include "Orbit_Number_Input.hpp"


int main() {
//MIGHT NEED TO ADD AN INPUT FOR THE LOGGER OUTPUT
    Logger *logger = Logger::getInstance("/Users/jonathonhope/Desktop/Cpp_424/Test_Logger_Files/System-Info.log");
    logger->log(PROGRAM_START,"");
    int measurement_set = Get_Measurement_Input();
    IOD state1;
    state1.Find_State(measurement_set, logger);
    
    //[ Logger call for the state of the spacecraft
    std::ostringstream message_state;
    message_state << "Position: [" << state1.r2_vec[0] << '\t' << state1.r2_vec[1] << '\t' << state1.r2_vec[2] <<"]" << " km" << '\n' << "Velocity: [" << state1.v2_vec[0] << '\t' << state1.v2_vec[1] << '\t' << state1.v2_vec[2] << "]" << " km/s" << '\n';
    logger->log(SC_STATE,message_state.str());
    
    //logger.log(SC_STATE,message_state);
    //]
    
        std::string orbit_type =  Get_Orbit_Type();
        int orbit_number = Get_Orbit_Number_Input(orbit_type);

    
    
    
    
    Arrival_Orbit arrival_orbit(orbit_number, orbit_type, logger);
    arrival_orbit.Choosing_orbit_type(orbit_number, orbit_type, logger);
  
    
    Create_Transfers(state1, arrival_orbit,logger);
    logger->log(PROGRAM_END,"");

    
    return 0;
    
    
    
}
