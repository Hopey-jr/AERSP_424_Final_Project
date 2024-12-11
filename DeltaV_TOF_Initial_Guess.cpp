//
//  DeltaV_TOF_Initial_Guess.cpp
//  Final_Project
//
//  Created by jonathon hope on 12/9/24.
//

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
