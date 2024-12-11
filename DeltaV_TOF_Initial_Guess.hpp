//
//  DeltaV_TOF_Initial_Guess.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef DeltaV_TOF_Initial_Guess_hpp
#define DeltaV_TOF_Initial_Guess_hpp

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
//#include "EOM_backwards_function.hpp"
//#include "EOM_STM_backwards_function.hpp"
#include "PushBackStateAndTime.hpp"

typedef std::vector<double> state_type;

// Observer class to handle integration control
class combined_observer {
public:
    combined_observer(boost::numeric::odeint::max_step_checker &checker, push_back_state_and_time &state_time_observer);

    void operator()(const state_type &w, double t);

private:
    boost::numeric::odeint::max_step_checker &m_checker;
    push_back_state_and_time &m_state_time_observer;
};

// Function to calculate Delta-V and TOF initial guess
std::vector<double> Get_DeltaV_TOF_Initial_Guess( std::vector<double> r_rot, std::vector<double> IC0, double Delta_V_upper, double Delta_V_lower, double TOF_upper, double TOF_lower, int n, double LU, double mu);
#endif /* DeltaV_TOF_Initial_Guess_hpp */
