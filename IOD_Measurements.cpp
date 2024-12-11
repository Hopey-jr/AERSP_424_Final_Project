//
//  IOD_Measurements.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

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




