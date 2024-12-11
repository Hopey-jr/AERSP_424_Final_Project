//
//  IOD_Measurements.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#ifndef IOD_Measurements_hpp
#define IOD_Measurements_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

class IOD_Measurements {
private:
   
    //std::ifstream file;
    std::ifstream file;
    std::ifstream file2;
    std::ifstream file3;
    std::ifstream file4;
    // Helper function to read data from a file
   // void readDataFromFile(const std::string& filepath, int measurement_number, double& value, int limit = 1);

public:
    double time;
    double azimuth;
    double elevation;
    std::vector<double> r_site;
    // Constructor
    IOD_Measurements();

    // Destructor
    ~IOD_Measurements();

// double time, double azimuth, double elevation,  std::vector<double> r_site, std::ifstream& file, std::ifstream& file2, std::ifstream& file3, std::ifstream& file4
    void Getting_IOD_Measurements(int measurement_number);

    // Accessor methods
    double getTime() const { return time; }
    double getAzimuth() const { return azimuth; }
    double getElevation() const { return elevation; }
    const std::vector<double>& getR_Site() const { return r_site; }
};
#endif /* IOD_Measurements_hpp */
