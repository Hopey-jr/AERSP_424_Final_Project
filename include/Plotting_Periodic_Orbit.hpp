//
//  Plotting_Periodic_Orbit.hpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#ifndef Plotting_Periodic_Orbit_hpp
#define Plotting_Periodic_Orbit_hpp
#include "Arrival_Orbit.hpp"

typedef std::vector< double > state_type;

std::vector< state_type > Plot_Periodic_Orbit(Arrival_Orbit& arrival_orbit);

#endif /* Plotting_Periodic_Orbit_hpp */
