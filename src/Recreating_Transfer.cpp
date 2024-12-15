//
//  Recreating_Transfer.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#include "Recreating_Transfer.hpp"
#include "Final_Outputs.hpp"
#include "IOD.hpp"
#include "Arrival_Orbit.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <boost/numeric/odeint.hpp>
#include "PushBackStateAndTime.hpp"
#include "EOM_function.hpp"
typedef std::vector< double > state_type;

std::vector< state_type > Recreate_Transfer(IOD& state1, Arrival_Orbit& arrival_orbit, double LU, double TU, double mu, double x_departure, double y_departure, double z_departure, double xdot_departure, double ydot_departure, double zdot_departure, double T1, double delta_vx, double delta_vy, double delta_vz, double T2_elapsed){
    
    namespace ublas = boost::numeric::ublas;
    using namespace std;
    using namespace boost::numeric::odeint;
    double abs_err = 1.0e-13; //absolute tolerance
    double rel_err = 1.0e-13; //relative tolerance
    state_type w(6);
    vector<state_type> w_vec1;
    vector<state_type> w_vec2;
    vector<double> times1;
    vector<double> times2;
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    auto controlled_stepper =make_controlled(abs_err,rel_err, error_stepper_type());
    w = {x_departure, y_departure, z_departure, xdot_departure, ydot_departure, zdot_departure};
    T1 = T1*3600/TU;
    double dt1 = (T1)/999;
    double T2 = T2_elapsed*3600/TU - T1;
    std::vector<double> times_vector1(1000);
    for (int p = 0; p < 1000; ++p) {
        times_vector1[p] = p * dt1;
    }
    
    //[ Integration to find the STM at one full time period
    integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop, w,  times_vector1.begin(),times_vector1.end(), dt1, push_back_state_and_time(w_vec1, times1));
    //]
    w.clear();
    
    w = {w_vec1[999][0],w_vec1[999][1],w_vec1[999][2],w_vec1[999][3]+delta_vx,w_vec1[999][4]+delta_vy,w_vec1[999][5]+delta_vz};
    double dt2 = (T2 )/999;
    std::vector<double> times_vector2(1000);
    for (int p = 0; p < 1000; ++p) {
        times_vector2[p] = p * dt2;
    }
    integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop, w,  times_vector2.begin(),times_vector2.end(), dt2, push_back_state_and_time(w_vec2, times2));
    
    
    w_vec1.insert(w_vec1.end(),w_vec2.begin(),w_vec2.end());
    return w_vec1;
    
}

