//
//  Orbit_Number_Input.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/15/24.
//

#include "Orbit_Number_Input.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
int Get_Orbit_Number_Input(std::string Orbit_Type) {
    
    int upper_bound = 0;
    std::string error_statement = "";
    std::ostringstream error_message;
    
    if (Orbit_Type == "LYAPUNOV_L1"){
        upper_bound = 1000;
    }else if(Orbit_Type == "LYAPUNOV_L2"){
        upper_bound = 1019;
    }else if(Orbit_Type == "HALO_L1_NORTHERN"){
        upper_bound = 1030;
    }else if(Orbit_Type == "HALO_L2_NORTHERN"){
        upper_bound = 761;
    }else if(Orbit_Type == "HALO_L1_SOUTHERN"){
        upper_bound = 972;
    }else if(Orbit_Type == "HALO_L2_SOUTHERN"){
        upper_bound = 761;
    }else{
        std::cout << "Invalid input" << std::endl;
    }
    
    error_message<<"Invalid input. Please enter a valid number of 1-" << upper_bound << std::endl;
    
    
    
    int value = 0;
    while (true) {
        try {
            std::cout << "What measurement set would you like? Choose 1-" << upper_bound << std::endl;
            std::cin >> value;

            // Check for invalid input (non-integer)
            if (std::cin.fail()) {
                throw std::logic_error(error_message.str());
            }

            // Check if the input is out of range
            if (value < 1 || value > upper_bound) {
                throw std::out_of_range(error_message.str());
            }

            // Valid input
            return value;
        } catch (const std::out_of_range& e) {
            std::cerr << e.what() << std::endl;

            // Clear the error state and ignore invalid input
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } catch (const std::logic_error& e) {
            std::cerr << e.what() << std::endl;

            // Clear the error state and ignore invalid input
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
}
