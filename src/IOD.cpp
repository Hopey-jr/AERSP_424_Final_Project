//
//  IOD.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "IOD.hpp"
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include "IOD_Measurements.hpp"
#include "helper.h"
#include "logger_V2.hpp"


// Constructor
IOD::IOD()  {
    
}

IOD::~IOD(){

}

void IOD::Find_State(int measurement_set,Logger *logger){
    logger->log(FUNCTION_START,"The position and velocity of the spacecraft will begin to be found");
    IOD_Measurements Measurement1;
    IOD_Measurements Measurement2;
    IOD_Measurements Measurement3;
    const int first_measurement = measurement_set*3-3;
    const int second_measurement =measurement_set*3-2;
    const int third_measurement = measurement_set*3-1;
    Measurement1.Getting_IOD_Measurements(first_measurement);
    Measurement2.Getting_IOD_Measurements(second_measurement);
    Measurement3.Getting_IOD_Measurements(third_measurement);
    
    //[ Logger call for the measurements
    std::ostringstream message_measurements;
    
    message_measurements << std::fixed << std::setprecision(6);
    message_measurements << '\t' << '\t' << '\t' << "Azimuth" << '\t' << "Elevation" << '\t' << "x_site" << '\t' << "y_site" << '\t' << "z_site" << '\t' << "time" << '\n' << "Measurement 1: " << Measurement1.azimuth << '\t' << Measurement1.elevation << '\t' << Measurement1.r_site[0] << '\t' << Measurement1.r_site[1] << '\t' << Measurement1.r_site[2] << '\t' << Measurement1.time << '\n' << "Measurement 2: " << Measurement2.azimuth << '\t' << Measurement2.elevation << '\t' << Measurement2.r_site[0] << '\t' << Measurement2.r_site[1] << '\t' << Measurement2.r_site[2] << '\t' << Measurement2.time << '\n' << "Measurement 3: " << Measurement3.azimuth << '\t' << Measurement3.elevation << '\t' << Measurement3.r_site[0] << '\t' << Measurement3.r_site[1] << '\t' << Measurement3.r_site[2] << '\t' << Measurement3.time << '\n';
     
    logger->log(SC_STATE,message_measurements.str());
    //]
    
    GaussOutput output;
    Gauss(Measurement1, Measurement2, Measurement3, mu, output);
    auto [R1_vec, R2_vec, R3_vec] = Get_position_vectors( output, mu);
    

     r2_vec = R2_vec;
     v2_vec = Gibbs( mu, R1_vec, R2_vec, R3_vec);
    logger->log(FUNCTION_END,"The position and velocity of the spacecraft have been to be found");

}




// Gauss method for initial orbit determination
std::vector<double> IOD::Gauss(const IOD_Measurements& Measurement1, const IOD_Measurements& Measurement2, const IOD_Measurements& Measurement3, double mu, GaussOutput& output) {
    double t1 = Measurement1.time;
    double t2 = Measurement2.time;
    double t3 = Measurement3.time;
    
    // Calculating L1, L2, L3
    Eigen::Vector3d L1(cos(Measurement1.azimuth) * cos(Measurement1.elevation), cos(Measurement1.elevation) * sin(Measurement1.azimuth), sin(Measurement1.elevation));
    Eigen::Vector3d L2(cos(Measurement2.azimuth) * cos(Measurement2.elevation), cos(Measurement2.elevation) * sin(Measurement2.azimuth), sin(Measurement2.elevation));
    Eigen::Vector3d L3(cos(Measurement3.azimuth) * cos(Measurement3.elevation), cos(Measurement3.elevation) * sin(Measurement3.azimuth), sin(Measurement3.elevation));
    
    // Calculating Ri
    Eigen::Vector3d R1(Measurement1.r_site[0], Measurement1.r_site[1], Measurement1.r_site[2]);
    Eigen::Vector3d R2(Measurement2.r_site[0], Measurement2.r_site[1], Measurement2.r_site[2]);
    Eigen::Vector3d R3(Measurement3.r_site[0], Measurement3.r_site[1], Measurement3.r_site[2]);
    
    // Calculating time intervals
    double tau1 = t1 - t2;
    double tau3 = t3 - t2;
    double tau13 = t3 - t1;
    
    // Triple Products
    double D0 = L1.dot(L2.cross(L3));
    
    double D11 = R1.dot(L2.cross(L3));
    double D21 = R2.dot(L2.cross(L3));
    double D31 = R3.dot(L2.cross(L3));
    
    double D12 = R1.dot(L1.cross(L3));
    double D22 = R2.dot(L1.cross(L3));
    double D32 = R3.dot(L1.cross(L3));
    
    double D13 = R1.dot(L1.cross(L2));
    double D23 = R2.dot(L1.cross(L2));
    double D33 = R3.dot(L1.cross(L2));
    
    // A and B coefficients
    double A = 1 / D0 * (-tau3 * D12 / tau13 + D22 + tau1 * D32 / tau13);
    double B = 1 / (6 * D0) * (-(pow(tau13, 2.0) - pow(tau3, 2.0)) * tau3 * D12 / tau13 + (pow(tau13, 2.0) - pow(tau1, 2.0)) * tau1 * D32 / tau13);
    
    // a, b, and c coefficients
    double a = -pow(A, 2.0) - 2 * A * L2.dot(R2) - pow(R2.norm(), 2.0);
    double b = -2 * mu * B * (A + L2.dot(R2));
    double c = -pow(mu, 2.0) * pow(B, 2.0);
    
    // Find the roots
    
    double r2= root_finder(a, b, c);
    output.r2  = r2;
    output.A = A;
    output.B = B;
    output.tau1 = tau1;
    output.tau3 = tau3;
    output.tau13 = tau13;
    output.D = { {0, D0, 0},{D11, D21, D31}, {D12, D22, D32}, {D13, D23, D33} };
    output.R = { {R1[0], R1[1], R1[2]}, {R2[0], R2[1], R2[2]}, {R3[0], R3[1], R3[2]} };
    output.L = { {L1[0], L1[1], L1[2]}, {L2[0], L2[1], L2[2]}, {L3[0], L3[1], L3[2]} };
    
    return {};
}

// Get position vectors
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> IOD::Get_position_vectors(const GaussOutput& gaussOutput, double mu) {
    
    double r2 = gaussOutput.r2;
        //double D0 = gaussOutput.D0;
        const auto& D = gaussOutput.D;
        const auto& R = gaussOutput.R;
        const auto& L = gaussOutput.L;
        double A = gaussOutput.A;
        double B = gaussOutput.B;
        double tau1 = gaussOutput.tau1;
        double tau3 = gaussOutput.tau3;
        double tau13 = gaussOutput.tau13;
    
    double p1 = 1 / D[0][1] * ((6 * (D[1][2] * tau1 / tau3 + D[1][1] * tau13 / tau3) * pow(r2, 3.0) + mu * D[1][2] * (pow(tau13, 2.0) - pow(tau1, 2.0)) * tau1 / tau3) / (6 * pow(r2, 3.0) + mu * (pow(tau13, 2.0) - pow(tau3, 2.0))) - D[1][0]);
    double p2 = A + mu * B * pow(r2, -3.0);
    double p3 = 1 / D[0][1] * ((6 * (D[3][0] * tau3 / tau1 - D[3][1] * tau13 / tau1) * pow(r2, 3.0) + mu * D[3][0] * (pow(tau13, 2.0) - pow(tau3, 2.0)) * tau3 / tau1) / (6 * pow(r2, 3.0) + mu * (pow(tau13, 2.0) - pow(tau1, 2.0))) - D[3][2]);
    
    std::vector<double> R1 = {R[0][0], R[0][1], R[0][2]};
    std::vector<double> R2 = {R[1][0], R[1][1], R[1][2]};
    std::vector<double> R3 = {R[2][0], R[2][1], R[2][2]};

    std::vector<double> L1 = {L[0][0], L[0][1], L[0][2]};
    std::vector<double> L2 = {L[1][0], L[1][1], L[1][2]};
    std::vector<double> L3 = {L[2][0], L[2][1], L[2][2]};
    
    
    std::vector<double> r1_vec(3);
    std::vector<double> r2_vec(3);
    std::vector<double> r3_vec(3);
    
    
    for(int i = 0; i < 3; ++i){
        r1_vec[i] = R1[i] + L1[i]*p1;
        r2_vec[i] = R2[i] + L2[i]*p2;
        r3_vec[i] = R3[i] + L3[i]*p3;
    }

    // Calculate position vectors (example for r1_vec, r2_vec)
    //std::vector<double> output = {};
    return std::make_tuple(r1_vec,r2_vec,r3_vec);
}

// Gibbs method for orbit determination
std::vector<double> IOD::Gibbs(double mu, const std::vector<double>& R1_vec, const std::vector<double>& R2_vec, const std::vector<double>& R3_vec) {
    Eigen::Vector3d r1_vec(R1_vec[0], R1_vec[1], R1_vec[2]);
    Eigen::Vector3d r2_vec(R2_vec[0], R2_vec[1], R2_vec[2]);
    Eigen::Vector3d r3_vec(R3_vec[0], R3_vec[1], R3_vec[2]);

    Eigen::Vector3d n = r1_vec.norm() * r2_vec.cross(r3_vec) + r2_vec.norm() * r3_vec.cross(r1_vec) + r3_vec.norm() * r1_vec.cross(r2_vec);
    Eigen::Vector3d d = r1_vec.cross(r2_vec) + r2_vec.cross(r3_vec) + r3_vec.cross(r1_vec);
    Eigen::Vector3d s = r1_vec * (r2_vec.norm() - r3_vec.norm()) + r2_vec * (r3_vec.norm() - r1_vec.norm()) + r3_vec * (r1_vec.norm() - r2_vec.norm());

    Eigen::Vector3d v = std::sqrt(mu / (n.norm() * d.norm())) * (d.cross(r2_vec) / r2_vec.norm() + s);
    std::vector<double> v2_vec = {v[0], v[1], v[2]};

    return v2_vec;
}

