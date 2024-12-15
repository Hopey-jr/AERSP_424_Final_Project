//Description:
//            The purpose of this class is to add attributes to an object, IOD_Measurements. The class will read from four text files to be able to add the attributes based on user input.

//Input:
//       [BASED ON USER INPUT] measurement_number - The number that correlates to the row of data the function needs to grab

//Output:
//       An object with the following attributes: azimuth, elevation, time, and r_site

//Local:
//       PI - Pi coming from the numebrs standard library
//       file(1-4) - The text file name that holds each attribute
//       term - The individual term from a single line
//       line - The string representation of the line chosen by the user for the function to read


//Function Call:
//           Measurement1.Getting_IOD_Measurements(first_measurement);


#include "IOD_Measurements.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <numbers>
#include <stdexcept>

// Constructor
IOD_Measurements::IOD_Measurements() {
}

// Destructor
IOD_Measurements::~IOD_Measurements() {
}

// Function to get IOD measurements
void IOD_Measurements::Getting_IOD_Measurements(int measurement_number) {
    double PI = std::numbers::pi;
    //[ Collecting the position of the sites
    file.open("/Users/jonathonhope/Desktop/Cpp 424/IOD_Initial_Conditions/R_SITES.txt");
    if (!file) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }
    std::string line;
    int i = 0;
    int row = 0;
    while (std::getline(file, line)) {
        if (row == measurement_number) {
            std::stringstream ss(line);
            double term;
            while (ss >> term && i < 3) {
                if (i < 3) {
                    r_site.push_back(term);
                }
                i = i + 1;
            }
            break;
        }
        row = row + 1;
    }
    file.close();
    line.clear();
    //]

    //[ Collecting the azimuth terms
    file2.open("/Users/jonathonhope/Desktop/Cpp 424/IOD_Initial_Conditions/AZIMUTH.txt");
    if (!file2) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    i = 0;
    row = 0;
    while (std::getline(file2, line)) {
        if (row == measurement_number) {
            std::stringstream ss(line);
            double term;
            while (ss >> term && i < 1) {
                azimuth = term*PI/180.0;
            }
        }
        row = row + 1;
    }
    file2.close();
    line.clear();
    //]

    //[ Collecting the elevation terms
    file3.open("/Users/jonathonhope/Desktop/Cpp 424/IOD_Initial_Conditions/ELEVATION.txt");
    if (!file3) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    i = 0;
    row = 0;
    while (std::getline(file3, line)) {
        if (row == measurement_number) {
            std::stringstream ss(line);
            double term;
            while (ss >> term && i < 1) {
                elevation = term*PI/180.0;
            }
        }
        row = row + 1;
    }
    file3.close();
    line.clear();
    //]

    //[ Collecting the time terms
    file4.open("/Users/jonathonhope/Desktop/Cpp 424/IOD_Initial_Conditions/TIMES.txt");
    if (!file4) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    i = 0;
    row = 0;
    while (std::getline(file4, line)) {
        if (row == measurement_number) {
            std::stringstream ss(line);
            double term;
            while (ss >> term && i < 1) {
                time = term;
            }
        }
        row = row + 1;
    }
    file4.close();
    //]
}




