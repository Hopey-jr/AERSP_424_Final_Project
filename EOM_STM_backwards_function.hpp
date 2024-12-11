//
//  EOM_STM_backwards_function.hpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

#ifndef EOM_STM_backwards_function_hpp
#define EOM_STM_backwards_function_hpp

#include <vector>

typedef std::vector<double> state_type;

void state_STM_prop_backwards(const state_type &w, state_type &dwdt, const double t);
#endif /* EOM_STM_backwards_function_hpp */
