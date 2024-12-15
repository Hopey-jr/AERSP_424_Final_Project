//
//  PushBackStateAndTime.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef PushBackStateAndTime_hpp
#define PushBackStateAndTime_hpp

#include <vector>

// Typedef for state_type
typedef std::vector<double> state_type;

// Struct to store states and times during integration
struct push_back_state_and_time {
    std::vector<state_type>& m_states;
    std::vector<double>& m_times;

    // Constructor
    push_back_state_and_time(std::vector<state_type>& states, std::vector<double>& times);

    // Overload operator() to push back state and time
    void operator()(const state_type& w, double t);
};
#endif /* PushBackStateAndTime_hpp */
