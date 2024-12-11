//
//  Arrival_Orbit.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef Arrival_Orbit_hpp
#define Arrival_Orbit_hpp


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

class Arrival_Orbit {
private:
    int orbit_number;
    std::string orbit_type;

public:
    int Arrival_orbit_number;
    std::string Orbit_Type;
    std::vector<double> Initial_Conditions;
    double Tp;
    std::ifstream file;

    Arrival_Orbit(int Arrival_orbit_number, std::string Orbit_Type);
    
    void Choosing_orbit_type(int Arrival_orbit_number, std::string Orbit_Type);
    void Getting_Initial_Conditions(int Arrival_orbit_number, std::ifstream& file);
};
#endif /* Arrival_Orbit_hpp */
