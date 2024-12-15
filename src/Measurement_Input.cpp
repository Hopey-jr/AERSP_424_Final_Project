//
//  Measurement_Input.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/15/24.
//

#include "Measurement_Input.hpp"

#include <iostream>
#include <stdexcept>

int Get_Measurement_Input() {
    int value = 0;
    while (true) {
        try {
            std::cout << "What measurement set would you like? Choose 1-8 " << std::endl;
            std::cin >> value;

            // Check for invalid input (non-integer)
            if (std::cin.fail()) {
                throw std::logic_error("Invalid input. Please enter a valid number of 1-8 '\n");
            }

            // Check if the input is out of range
            if (value < 1 || value > 8) {
                throw std::out_of_range("Invalid input. Please enter a valid number of 1-8 '\n");
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
