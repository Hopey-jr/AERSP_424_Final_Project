//  Description:
//               The purpose of this function is to determine which branches of the manifold intersect with the Moon. For the continuation process it is helpful to keep all the branches of the manifold so there are no big jumps of initial conditions. Once the continuation procoess is finished, all the branches that originally hit the manifold will be excluded from being saved to the text file.


//Inputs:
//       stable_manifold - The 2 dimensional vector that holds the initial conditions of the stable manifold, along with the time it takes to hit the manifold
//       hits_moon - The vector that holds which branches of the manifold hit the Moon
//       start - The value that indicates which chuck of the processor to start with
//       end - The value that indicated which chunk of the processor to end with
//       abs_tol - The absolute tolerance for the ode intergration
//       rel_tol - The relative tolerance for the ode integration
//       LU - The non-dimensional distance unit equal to distance from the Earth to the Moon
//       mu - The non-dimensional reduced mass of the Earth-Moon system


//Output:
//       hits_moon - The vector that holds which branches of the manifold hit the Moon


//Locals:
//       w - The vector and variable used to hold the initial conditions for the integration that include both the initial state and initial STM
//       dt - The initial timestep used for the adaptive timestep integration
//       w_vec - The state and value of the STM from 0 to tau
//       times - The timestamps of each integration
//       times_vector - The vector that holds the time span for the ode function
//       TOF - The time of flight along the manifold
//       distance_moon - The distance to the Moon in non-dimensional units
//       local_hits - A counter that keeps a tally of how many manifolds hit the Moon

//Dependent Functions:
//       PushBackStateAndTime.hpp" - The purpose of this function is to save the time and state in their respective vectors during the ode integration
//       EOM_function.hpp" - The purpose of this function is to propagate the equations of motion of the CR3BP in the ode function

//Function call: 
//       Check_If_Hits_Moon( std::vector<std::vector<double>>& stable_manifold, Check_If_Hits_Moon( hits_moon,  start,  end,  abs_err,  rel_err,  LU,  mu);




#include "Hits_Moon_Check.hpp"
#include <thread>
#include <vector>
#include <mutex>
#include <cmath>
#include <algorithm>
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
void Check_If_Hits_Moon( std::vector<std::vector<double>>& stable_manifold,
                         std::vector<double>& hits_moon,
                   int start, int end,
                   double abs_err, double rel_err, double LU, double mu) {

    int local_hits = 0;
    std::vector<double> w(6);
    std::vector<state_type> w_vec;
    std::vector<double> times;

    
    for (int j = start; j < end; ++j) {
        
        
        
        
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

       //thread_local std::vector<std::vector<double>> w_vec;
       //thread_local std::vector<double> times;

        // Perform integration
        integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()),
                        state_prop, w, times_vector.begin(), times_vector.end(), dt,
                        push_back_state_and_time(w_vec, times));

        for (int moon_test = 0; moon_test < w_vec.size(); ++moon_test) {
            double distance_moon = std::pow(std::pow((w_vec[moon_test][0] - (1.0 - mu)), 2.0) +
                                            std::pow(w_vec[moon_test][1], 2.0) +
                                            std::pow(w_vec[moon_test][2], 2.0), 0.5);

            if (distance_moon < 1740 / LU) {
                
                    std::lock_guard<std::mutex> lock(hits_moon_mutex);
                    hits_moon.push_back(j);
                
                ++local_hits;
                break;
            }
        }
        
        
        
        
    }


}
