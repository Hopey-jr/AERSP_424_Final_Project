//
//  Arrival_Orbit.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "Arrival_Orbit.hpp"
Arrival_Orbit::Arrival_Orbit(int Arrival_orbit_number, std::string Orbit_Type)
    : orbit_number(Arrival_orbit_number), orbit_type(Orbit_Type) {
    this->Arrival_orbit_number = Arrival_orbit_number;
    this->Orbit_Type = Orbit_Type;
    std::cout << "Arrival Orbit created of type " << Orbit_Type << " with orbit number " << Arrival_orbit_number << std::endl;
}

void Arrival_Orbit::Choosing_orbit_type(int Arrival_orbit_number, std::string Orbit_Type) {
    if (Orbit_Type == "LYAPUNOV_L1") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Lyapunov_L1");
    } else if (Orbit_Type == "LYAPUNOV_L2") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Lyapunov_L2");
    } else if (Orbit_Type == "HALO_L1_NORTHERN") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Halo_L1_Northern");
    } else if (Orbit_Type == "HALO_L2_NORTHERN") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Halo_L2_Northern");
    } else if (Orbit_Type == "HALO_L1_SOUTHERN") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Halo_L1_Southern");
    } else if (Orbit_Type == "HALO_L2_SOUTHERN") {
        file.open("/Users/jonathonhope/Desktop/Graduate_Research/Complete_Catelogs/Halo_L2_Southern");
    } else {
        std::cerr << "Invalid orbit type" << std::endl;
        return;
    }

    if (!file) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    Getting_Initial_Conditions(Arrival_orbit_number, file);
    file.close();
}

void Arrival_Orbit::Getting_Initial_Conditions(int Arrival_orbit_number, std::ifstream& file) {
    std::string line;
    int row = 0;
    int i = 0;

    while (std::getline(file, line)) {
        row = row + 1;
        if (row == Arrival_orbit_number) {
            std::stringstream ss(line);
            double term;
            while (ss >> term && i < 16) {
                if (i < 6) {
                    Initial_Conditions.push_back(term);
                }
                if (i == 12) {
                    Tp = term;
                }
                i = i + 1;
            }
            break;
        }
    }
}
