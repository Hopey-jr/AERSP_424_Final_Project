//  Description:
//               The purpose of this code is to find the closest point to the departure orbit given an array of TOF and delta_v guesses

//Inputs:
//       r_rot - The position of the satellite that will be used as the constraint in finding the optimal delta_v and Tp
//       IC0 - The initial conditions for the branch of the manifold coming from the 'manifold_positive' or 'manifold_negative' functions
//       Delta_V_upper - The upper bound of the delta_v guess
//       Delta_V_lower - The lower bound of the delta_v guess
//       TOF_upper - The upper bound of the time of flight guess
//       TOF_lower - The lower bound of the time of flight guess
//       n - the number of values you want to guess for for delta_v and TOF upper and lower
//       LU - The characteristic length unit derived from the average distance between the Earth and Moon
//       mu - The nondimensional reduced mass of the Earth-Moon system [LU^3/TU^2]
//Output:
//       output [3x1]:
//                      Tp - The time that corresponds to when the trajectroy hits poincare section
//                      delta_V - the delta_v to when the trajectroy hits poincare section
//                      minval - The value of the minimum error

//Locals:

//       abs_err - The absolute tolerance for the error
//       rel_err - The relative tolerance for the error
//       w_vec - The state and value of the STM from t0 to tf
//       times - The timestamps of each integration
//       dt_deltaV - The stepsize of delta_V between Delta_V_lower and Delta_V_upper
//       dt_Tps - The stepsize of TOF between TOF_lower and TOF_upper
//       delta_v - The current value of delta_v for testing
//       TOF_Guess - The current value of TOF for testing
//       r - The initial position passed into the function
//       v - The initial velocity passed into the function
//       v_unit_vector - The direction of the velocity vector
//       deltaV_vector - The vector of all the delta_V values used for guessing
//       Tps_vector - The vector of all the TOF values used for guessing
//       v_unit_vector - The direction of the initial velocity vector
//       v_transfer - The velocity after the delta_v is applied
//       w - The vector and variable used to hold the initital conditions for the integration that include both the initital state and initial STM
//       dt_Guess - The timestep used for the integration for each guess of delta_v and TOF
//       times_vector_Guess - The time vector used for each guess of delta_v and TOF
//       end_of_int2 - The integer value that holds the last position in v_vec and times
//       te_hit - The vector that holds the timestamp if the trajectory hits the poincare section
//       delta_v_hit - The vector that holds the corresponding deltaV if the trajectory hits the poincare section
//       q - Counter for the amount of trajectories that hit the
//       t_end - The final time it is allowed to integrate to. If the integration hits this point, the trajectory does not hit the poincare section
//       dt - The change in time for the integration
//       step_count - The counter for the amount of steps each integration is taking to converge
//       max_steps - The amount of steps the integration is allowed to take before it stops
//       edge_tol - The allowable tolerance that the trajectory must be considered close to the arrival orbit
//       min - The position in the te_hit and delta_v_hit vectors that describes the minimum value in te_hit vector
//       min_val - The minimum value of the te_hit vector

//Dependent Functions:
//       EOM_backwards_function.h - Propagates the STM and state from tf to t0
//       PushBackStateAndTime.h - Creates the arrays to hold the state and time

//function call:
//     std::vector<double> Output = Get_DeltaV_TOF_IC_V2(r_rot, IC0,  delta_V_high,  delta_V_low,  TOF_high,  TOF_low, n, LU, mu);

#include "DeltaV_TOF_Initial_Guess.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include "EOM_backwards_function.hpp"
#include "EOM_STM_backwards_function.hpp"
#include "PushBackStateAndTime.hpp"

typedef std::vector<double> state_type;
namespace ublas = boost::numeric::ublas;

combined_observer::combined_observer(boost::numeric::odeint::max_step_checker &checker, push_back_state_and_time &state_time_observer)
    : m_checker(checker), m_state_time_observer(state_time_observer) {}

void combined_observer::operator()(const state_type &w, double t) {
    m_checker();
    m_state_time_observer(w, t);
}

std::vector<double> Get_DeltaV_TOF_Initial_Guess(std::vector<double> r_rot, std::vector<double> IC0, double Delta_V_upper, double Delta_V_lower,double TOF_upper, double TOF_lower,int n, double LU, double mu) {
    using namespace boost::numeric::odeint;

    state_type w(6);
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    auto controlled_stepper = make_controlled(1e-13, 1e-13, error_stepper_type());

    std::vector<state_type> w_vec;
    std::vector<double> times;

    double dt_deltaV = (Delta_V_upper - Delta_V_lower) / n;
    double dt_Tps = (TOF_upper - TOF_lower) / n;

    std::vector<double> deltaV_vector(n), Tps_vector(n);
    Eigen::MatrixXd min_value(n, n);

    for (int m = 0; m < n; ++m) {
        deltaV_vector[m] = m * dt_deltaV + Delta_V_lower;
        Tps_vector[m] = m * dt_Tps + TOF_lower;
    }

    Eigen::Vector3d r(IC0[0], IC0[1], IC0[2]);
    Eigen::Vector3d v(IC0[3], IC0[4], IC0[5]);
    Eigen::Vector3d v_unit_vector = v.normalized();

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            std::fill(w.begin(), w.end(), 0);
            w_vec.clear();
            times.clear();

            double delta_v = deltaV_vector[i];
            double TOF_Guess = Tps_vector[j];
            Eigen::Vector3d v_transfer = delta_v * v_unit_vector;

            w[0] = r[0];
            w[1] = r[1];
            w[2] = r[2];
            w[3] = v_transfer[0] + v[0];
            w[4] = v_transfer[1] + v[1];
            w[5] = v_transfer[2] + v[2];

            std::vector<double> target_times(1000);
            double t_start = 0.0, t_end = TOF_Guess;
            for (size_t k = 0; k < target_times.size(); ++k) {
                target_times[k] = t_start + k * (t_end - t_start) / (target_times.size() - 1);
            }

            double t = t_start;
            size_t step_count = 0, max_steps = 1000, next_time_idx = 0;
            while (t < t_end && step_count < max_steps && next_time_idx < target_times.size()) {
                double dt = target_times[next_time_idx] - t;

                controlled_step_result result = controlled_stepper.try_step(state_prop_backwards, w, t, dt);
                if (result == success) {
                    if (std::abs(t - target_times[next_time_idx]) < 1e-12) {
                        w_vec.push_back(w);
                        times.push_back(t);
                        next_time_idx++;
                    }
                }
                step_count++;
            }

            int end_of_int = w_vec.size() - 1;
            double distance = 0.0;

            
                distance = std::sqrt(std::pow(w_vec[end_of_int][0] - r_rot[0], 2.0) +
                                     std::pow(w_vec[end_of_int][1]-r_rot[1], 2.0) +
                                     std::pow(w_vec[end_of_int][2]-r_rot[2], 2.0));

            min_value(j, i) = std::abs(distance);
        }
    }

    Eigen::Index Row, Col;
    double minval = min_value.minCoeff(&Row, &Col);

    std::vector<double> output = { Tps_vector[Row], deltaV_vector[Col], minval };
    return output;
}
