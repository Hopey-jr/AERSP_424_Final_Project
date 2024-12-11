//
//  EOM_function.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#ifndef EOM_function_hpp
#define EOM_function_hpp

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

// Typedef for state_type
typedef std::vector<double> state_type;

// Namespace alias for Boost uBLAS
namespace ublas = boost::numeric::ublas;

// Function declaration
void state_prop(const state_type &w, state_type &dwdt, const double t);

#endif /* EOM_function_hpp */
