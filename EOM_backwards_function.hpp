//
//  EOM_backwards_function.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef EOM_backwards_function_hpp
#define EOM_backwards_function_hpp
#include <vector> 
#include <boost/numeric/odeint.hpp>

typedef std::vector<double> state_type;

void state_prop_backwards(const state_type &w, state_type &dwdt, const double t);
#endif /* EOM_backwards_function_hpp */
