//
//  main.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/9/24.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <chrono>
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

int main() {
    std::cout << "What measurement set would you like? Choose 1-8" << std::endl;
    int measurement_set;
    std::cin >> measurement_set;
    IOD state1;
   state1.Find_State(measurement_set);
    std::cout << "Based on the chosen measurement set, the position and velocity of the spacecraft are: " << std::endl;
    std::cout << "["<< state1.r2_vec[0] << '\t' << state1.r2_vec[1] << '\t' << state1.r2_vec[2] <<"]" << " km" <<std::endl;
    std::cout << "["<< state1.v2_vec[0] << '\t' << state1.v2_vec[1] << '\t' << state1.v2_vec[2] << "]" << " km/s"<< std::endl;
    
    
    
    
    std::cout << "What arrival type would you like?" << std::endl;
    std::string orbit_type;
    std::cin >> orbit_type;
    
    
    std::cout << "What arrival orbit number would you like?" << std::endl;
    int orbit_number;
    std::cin >> orbit_number;
    
    Arrival_Orbit arrival_orbit(orbit_number, orbit_type);
    arrival_orbit.Choosing_orbit_type(orbit_number, orbit_type);
  
    
    Create_Transfers(state1, arrival_orbit);
    
    
    return 0;
    
    
    
}

// L1 Lyap, orbit 40. Switch to 40 iterations max. also switch up the outputs
