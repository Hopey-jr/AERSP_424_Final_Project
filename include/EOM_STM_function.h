

#ifndef EOM_STM_function_h
#define EOM_STM_function_h
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

// Typedef for state_type
typedef std::vector<double> state_type;

// Namespace alias for Boost uBLAS
namespace ublas = boost::numeric::ublas;

// Function declaration
void state_STM_prop(const state_type &w, state_type &dwdt, const double t);

#endif /* EOM_STM_function_h */
