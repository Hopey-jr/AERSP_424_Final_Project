//  Description:
//               The purpose of this class is to take the user input of orbit type and orbit number, read from catalog text files, and save the attributes of the initial conditions and orbital period to an object


//Inputs:
//       Arrival_orbit_number - The number orbit the user wants from the catalog
//       Orbit_Type - The orbit family the user chooses from the range of catalogs

//Output:
//       An object with the following attributes:
//          Initial_Conditions - The state (position and velocity) for the initial conditions to produce the orbit in non-dimensional units
//          Tp - The orbital period of the orbit in non-dimensional units
//Locals:
//       file - The file that holds all the initial conditions depending on the user input

//Dependent Functions:
//

//Function call:


#include "Arrival_Orbit.hpp"

Arrival_Orbit::Arrival_Orbit(int Arrival_orbit_number, std::string Orbit_Type, Logger *logger)
    : orbit_number(Arrival_orbit_number), orbit_type(Orbit_Type) {
    this->Arrival_orbit_number = Arrival_orbit_number;
    this->Orbit_Type = Orbit_Type;
    std::cout << "Arrival Orbit created of type " << Orbit_Type << " with orbit number " << Arrival_orbit_number << std::endl;
}

void Arrival_Orbit::Choosing_orbit_type(int Arrival_orbit_number, std::string Orbit_Type, Logger *logger) {
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

    Getting_Initial_Conditions(Arrival_orbit_number, file, logger);
    file.close();
}

void Arrival_Orbit::Getting_Initial_Conditions(int Arrival_orbit_number, std::ifstream& file, Logger *logger) {
    std::string line;
    int row = 0;
    int i = 0;
    double C = 0;
    double stab = 0;
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
                if (i == 14){
                     C = term;
                }
                if (i == 15){
                     stab = term;
                }
                i = i + 1;
            }
            break;
        }
    }
    std::ostringstream message_PO;
    message_PO << "The initial conditions of the orbit is: [" << Initial_Conditions[0] << ", " << Initial_Conditions[1] << ", " << Initial_Conditions[2] << ", " << Initial_Conditions[3] << ", " << Initial_Conditions[4] << ", " << Initial_Conditions[5] << "] [LU | LU/TU]" << '\n' << "The orbital period is: " << Tp << "TU" << '\n' << "The Jacobi energy is: " << C << '\n' << "The stability is: " << stab << '\n';
    logger->log(SC_STATE,message_PO.str());
}
