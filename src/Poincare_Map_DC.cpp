//  Description:
//               The purpose of this code is to find the state when hitting the poincare section using differential corrections. This code will have a similar result to the 'events' function in matlab

//Inputs:
//       ICs - The initial conditions for the branch of the manifold coming from the 'manifold_positive' or 'manifold_negative' functions
//       poincare_section_x - The final position for which the differential corrector is shooting for
//       mu - The nondimensional reduced mass of the Earth-Moon system [LU^3/TU^2]
//Output:
//       output [7]:
//                      State - The state when the trajectory hits the poincare section
//                      TOF - The time spent along the manifold until it hits the poincare section


//Locals:
//       tau - The time used to find the initial good guess for the corrector to ensure the trajectory passes the poincare section
//       closest_time - The timestamp that is is the closest approach for a single trajectroy to find a good guess for the differential corrector
//       closest_appraoch - The corresponding position to see how close the closest approach is to the given poincare section
//       err - The difference between the poincare section and the x-position at the end of integration
//       tol - The tolerance the differential corrector must hit for the solution to be considered converged. For this case, in the interest of computational power and time, the tol is set to ~10^-9
//       w - The vector and variable used to hold the initital conditions for the integration that include both the initital state and initial STM
//       w_vec - The state and value of the STM from t0 to tf
//       times - The timestamps of each integration
//       x_final - The x-position at the end of the integration
//       xdot_final - The x-component of velocity at the end of the integration
//       constraint - The constraint based on the poincare section
//       DF - The correction factor in the differential corrector



//Dependent Functions:
//       EOM_STM_backwards_function.h - Propagates the STM and state from tf to t0
//       EOM_backwards_function.h - - Propagates the state from tf to t0
//       PushBackStateAndTime.h - Creates the arrays to hold the state and time

//function call:
//             std::vector<double> Output = Get_Poincare_Section_2( mu,ICs,poincare_section_x);

#include "Poincare_Map_DC.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <complex>
#include <numbers>
#include <utility>
#include "EOM_STM_backwards_function.hpp"
#include "EOM_backwards_function.hpp"
#include "PushBackStateAndTime.hpp"

typedef std::vector<double> state_type;
namespace ublas = boost::numeric::ublas;

std::vector<double> Get_Poincare_Section_2(double mu, std::vector<double> ICs, double poincare_section_x) {
    using namespace std;
    using namespace boost::numeric::odeint;

    double abs_err = 1.0e-13; // Absolute tolerance
    double rel_err = 1.0e-13; // Relative tolerance

    vector<state_type> w_vec;
    vector<double> times;

    state_type w(6);
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;

    double tau = 12.0;
    double dt = tau / 3999.0;

    std::fill(w.begin(), w.end(), 0);

    w[0] = ICs[0];
    w[1] = ICs[1];
    w[2] = ICs[2];
    w[3] = ICs[3];
    w[4] = ICs[4];
    w[5] = ICs[5];

    w_vec.clear();
    times.clear();

    integrate_adaptive(
        make_controlled(abs_err, rel_err, error_stepper_type()),
        state_prop_backwards,
        w,
        0.0,
        tau,
        dt,
        push_back_state_and_time(w_vec, times)
    );

    double closest_time = 0;
    double closest_approach = 100;
    double time_tol = pow(10, -3);
    for (int k = 0; k < w_vec.size(); k++) {
        if (std::abs(w_vec[k][0] - poincare_section_x) < time_tol) {
            closest_time = times[k];
            closest_approach = std::abs(w_vec[k][0] - poincare_section_x);
            break;
        }
    }

    double TOF = closest_time;
    double err = 1;
    double tol = pow(10, -7);
    std::vector<double> Output(7);
    std::fill(Output.begin(), Output.end(), 0);

    for (int h = 0; h < 20; ++h) {
        std::fill(w.begin(), w.end(), 0);
        w_vec.clear();
        double dt2 = TOF / 1999.0;
        std::vector<double> times_vector2(2000);

        for (size_t m = 0; m < times_vector2.size(); ++m) {
            times_vector2[m] = m * dt2;
        }

        w[0] = ICs[0];
        w[1] = ICs[1];
        w[2] = ICs[2];
        w[3] = ICs[3];
        w[4] = ICs[4];
        w[5] = ICs[5];

        times.clear();
        integrate_times(
            make_controlled(abs_err, rel_err, error_stepper_type()),
            state_prop_backwards,
            w,
            times_vector2.begin(),
            times_vector2.end(),
            dt2,
            push_back_state_and_time(w_vec, times)
        );

        int end_of_int2 = w_vec.size() - 1;
        double x_final = w_vec[end_of_int2][0];
        double xdot_final = w_vec[end_of_int2][3];

        double constraint = x_final - poincare_section_x;
        double DF = constraint / xdot_final;
        err = std::abs(constraint);

        if (err < tol) {
            Output[0] = w_vec[end_of_int2][0];
            Output[1] = w_vec[end_of_int2][1];
            Output[2] = w_vec[end_of_int2][2];
            Output[3] = w_vec[end_of_int2][3];
            Output[4] = w_vec[end_of_int2][4];
            Output[5] = w_vec[end_of_int2][5];
            Output[6] = TOF;
            return Output;
        }

        TOF = TOF + 0.9*DF;
        if (TOF > 12) {
            break;
        }
    }

    return Output;
}
