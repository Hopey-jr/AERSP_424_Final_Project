//
//  Hits_Moon_Check.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/11/24.
//

#include "Hits_Moon_Check.hpp"
#include <thread>
#include <vector>
#include <mutex>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "PushBackStateAndTime.hpp"
#include "EOM_function.hpp"

// Mutex for shared resources
std::mutex hits_moon_mutex;
using namespace boost::numeric::odeint;
typedef std::vector<double> state_type;
typedef runge_kutta_cash_karp54<std::vector<double>> error_stepper_type;

// Function to process a chunk of the stable manifold
void Check_If_Hits_Moon(const std::vector<std::vector<double>>& stable_manifold,
                   std::vector<int>& hits_moon,
                   int& number_of_hits,
                   size_t start, size_t end,
                   double abs_err, double rel_err, double LU, double mu) {

    int local_hits = 0;
    std::vector<double> w(6);
    for (size_t j = start; j < end; ++j) {
        std::fill(w.begin(), w.end(), 0);
        for (size_t i = 0; i < 6; ++i) {
            w[i] = stable_manifold[j][i];
        }
        double TOF = stable_manifold[j][6];
        double dt = TOF / 4999.0;
        std::vector<double> times_vector(5000);
        for (size_t p = 0; p < 5000; ++p) {
            times_vector[p] = p * dt;
        }

       thread_local std::vector<std::vector<double>> w_vec;
       thread_local std::vector<double> times;

        // Perform integration
        integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()),
                        state_prop, w, times_vector.begin(), times_vector.end(), dt,
                        push_back_state_and_time(w_vec, times));

        for (size_t moon_test = 0; moon_test < w_vec.size(); ++moon_test) {
            double distance_moon = std::pow(std::pow((w_vec[moon_test][0] - (1.0 - mu)), 2.0) +
                                            std::pow(w_vec[moon_test][1], 2.0) +
                                            std::pow(w_vec[moon_test][2], 2.0), 0.5);

            if (distance_moon < 1740 / LU) {
                {
                    std::lock_guard<std::mutex> lock(hits_moon_mutex);
                    hits_moon.push_back(j);
                }
                ++local_hits;
                break;
            }
        }
    }

    // Update global number of hits
    {
        std::lock_guard<std::mutex> lock(hits_moon_mutex);
        number_of_hits += local_hits;
    }
}
