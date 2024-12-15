//
//  Orbit_Type_Input.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/15/24.
//

#include "Orbit_Type_Input.hpp"

#include <iostream>
#include <set>
#include <string>
#include <stdexcept>

std::string Get_Orbit_Type() {
    // Define the valid strings
    const std::set<std::string> validInputs = {
        "LYAPUNOV_L1",
        "LYAPUNOV_L2",
        "HALO_L1_NORTHERN",
        "HALO_L2_NORTHERN",
        "HALO_L1_SOUTHERN",
        "HALO_L2_SOUTHERN"
    };

    std::string input;
    while (true) {
        try {
            std::cout << "What arrival type would you like? Enter one of the following options: LYAPUNOV_L1, LYAPUNOV_L2, HALO_L1_NORTHERN, HALO_L2_NORTHERN, HALO_L1_SOUTHERN, HALO_L2_SOUTHERN" << std::endl;
            std::cin >> input;

            // Check if input is in the validInputs set
            if (validInputs.find(input) == validInputs.end()) {
                throw std::out_of_range("Invalid input. Please try again from these options: LYAPUNOV_L1, LYAPUNOV_L2, HALO_L1_NORTHERN, HALO_L2_NORTHERN, HALO_L1_SOUTHERN, HALO_L2_SOUTHERN '\n' ");
            }

            // Valid input
            return input;
        } catch (const std::out_of_range& e) {
            std::cerr << e.what() << std::endl;

            // Clear the input buffer to prevent infinite loops
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
}
