//
//  PushBackStateAndTime.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#include "PushBackStateAndTime.hpp"
// Constructor implementation
push_back_state_and_time::push_back_state_and_time(std::vector<state_type>& states, std::vector<double>& times)
    : m_states(states), m_times(times) {}

// Operator() implementation
void push_back_state_and_time::operator()(const state_type& w, double t) {
    m_states.push_back(w);
    m_times.push_back(t);
}
