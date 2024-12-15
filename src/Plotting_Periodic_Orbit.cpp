//
//  Plotting_Periodic_Orbit.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#include "Plotting_Periodic_Orbit.hpp"
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

std::vector< state_type > Plot_Periodic_Orbit(Arrival_Orbit& arrival_orbit){
    
    namespace ublas = boost::numeric::ublas;
    using namespace std;
    using namespace boost::numeric::odeint;
    double abs_err = 1.0e-13; //absolute tolerance
    double rel_err = 1.0e-13; //relative tolerance
    state_type w(6);
    vector<state_type> w_vec1;
    vector<double> times1;

    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    auto controlled_stepper =make_controlled(abs_err,rel_err, error_stepper_type());
    w = {arrival_orbit.Initial_Conditions[0], arrival_orbit.Initial_Conditions[1], arrival_orbit.Initial_Conditions[2], arrival_orbit.Initial_Conditions[3], arrival_orbit.Initial_Conditions[4], arrival_orbit.Initial_Conditions[5]};
    double dt1 = arrival_orbit.Tp/999;
    std::vector<double> times_vector1(1000);
    for (int p = 0; p < 1000; ++p) {
        times_vector1[p] = p * dt1;
    }
    
    //[ Integration to find the STM at one full time period
    integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop, w,  times_vector1.begin(),times_vector1.end(), dt1, push_back_state_and_time(w_vec1, times1));
    
    return w_vec1;
}
