//
//  IOD.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef IOD_hpp
#define IOD_hpp


#include <vector>
#include <Eigen/Dense>
#include <tuple>
#include "RootFinder.hpp"
#include "IOD_Measurements.hpp"
#include "helper.h"
#include "logger_V2.hpp"
struct GaussOutput {
    double D0;
    std::vector<std::vector<double>> D;
    std::vector<std::vector<double>> R;
    std::vector<std::vector<double>> L;
    double A;
    double B;
    double tau1;
    double tau3;
    double tau13;
    double r2;
};



// Class for Initial Orbit Determination (IOD) methods
class IOD {
    
    
private:
   const double mu = 3.986*pow(10,5);  // Gravitational parameter (mu)
    
public:
    std::vector<double> r2_vec;
    std::vector<double> v2_vec;  
    
    // Constructor
    IOD();

    //Destructor
    ~IOD();
    
   void Find_State(int measurement_set,Logger *logger);
    
    // Method to find roots using Gauss' method
    std::vector<double> Gauss(const IOD_Measurements& Measurement1, const IOD_Measurements& Measurement2, const IOD_Measurements& Measurement3, double mu, GaussOutput& output);

    // Method to get position vectors
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Get_position_vectors(const GaussOutput& output, double mu);

    // Method to apply Gibbs' method for orbit determination
    std::vector<double> Gibbs(double mu, const std::vector<double>& R1_vec, const std::vector<double>& R2_vec, const std::vector<double>& R3_vec);

    // Getter for the position and velocity vectors
   // const std::vector<double>& getPositionVector()
   // const std::vector<double>& getVelocityVector()




};

#endif /* IOD_hpp */
