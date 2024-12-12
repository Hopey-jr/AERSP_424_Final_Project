//
//  Transfers.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#include "Transfers.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <complex>
#include <numbers>
#include <utility>
#include <thread>

//#include "Hits_Moon_Check.hpp"
#include "DeltaV_TOF_Initial_Guess.hpp"
#include "Manifold_Negative.hpp"
#include "Manifold_Positive.hpp"
#include "Poincare_Map_DC.hpp"
#include "EOM_backwards_function.hpp"
#include "EOM_STM_backwards_function.hpp"
#include "EOM_STM_function.h"
#include "PushBackStateAndTime.hpp"
#include "EOM_function.hpp"
#include "ECI2Rotating_Frame_Coord_Transformation.hpp"
#include "PseudoInverse.hpp"
#include "Final_Outputs.hpp"
//#include "Proccessing_Negative_Manifold.hpp"
//#include "Proccessing_Positive_Manifold.hpp"



void Create_Transfers(IOD& state1, Arrival_Orbit& arrival_orbit){

    //[ Definitions and constants
    typedef std::vector< double > state_type;
    namespace ublas = boost::numeric::ublas;
     const double mu = 1.215058560962404E-2;
    using namespace std;
    using namespace boost::numeric::odeint;
    double abs_err = 1.0e-13; //absolute tolerance
    double rel_err = 1.0e-13; //relative tolerance
    state_type w(42);
    vector<state_type> w_vec;
    vector<state_type> w_vec2;
    vector<state_type> w_vec3;
    vector<double> times;
    double LU_JPL = 389703;
    double TU_JPL = 382981;
    double mu_E = pow(LU_JPL,3)/pow(TU_JPL,2);
    double LU = 384400;
    double TU = 3.751917661977041e+05;
    double hour = TU/3600.0;
    double j = 1;
    double Tp;
    double Tp_depart;
    double TOF;
    double delta_V;
    double v2_0;
    double Tps;
    double Jac2_dep;
    const int number_of_manifolds = 150;
    std::vector<double> delta_vx;
    std::vector<double> delta_vy;
    std::vector<double> delta_vz;
    std::vector<double> x_departure;
    std::vector<double> y_departure;
    std::vector<double> z_departure;
    std::vector<double> xdot_departure;
    std::vector<double> ydot_departure;
    std::vector<double> zdot_departure;
    std::vector<double> T1;
    std:vector<double> T2_elapsed;
    std::vector<double> T2;
    std::vector<double> hits_moon;
    int number_of_manifolds_at_poincare = 0;
    double tol = pow(10,-12);
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    //]
    
    
    //[ Finding the stable manifold
    std::vector<double> Initial_Conditions = arrival_orbit.Initial_Conditions;
    Tp = arrival_orbit.Tp;
    std::vector<double> ICs(6);
    
    
    /*
    std::vector<std::vector<double>> stable_negative_state(number_of_manifolds, std::vector<double>(7, 100));
        std::vector<std::vector<double>> stable_positive_state(number_of_manifolds, std::vector<double>(7, 100));
    
    
    auto arrays_neg = Generate_Negative_Manifolds( Initial_Conditions,  Tp,  mu,  number_of_manifolds);
    auto arrays_pos = Generate_Positive_Manifolds( Initial_Conditions,  Tp,  mu,  number_of_manifolds);
    std::vector<std::vector<double>> ICs_Stable_Negative = arrays_neg.first;
    std::vector<std::vector<double>> ICs_Unstable_Negative = arrays_neg.second;
    std::vector<std::vector<double>> ICs_Stable_Positive = arrays_pos.first;
    std::vector<std::vector<double>> ICs_Unstable_Positive = arrays_pos.second;
    double poincare_section_x = -mu;
        int hit_poincare_section = 0;

        // Split the work into multiple threads for negative and positive manifolds
        int num_threads = std::thread::hardware_concurrency(); // Use available hardware threads
        int chunk_size_neg = ICs_Stable_Negative.size() / num_threads;
        int chunk_size_pos = ICs_Stable_Positive.size() / num_threads;
        
        std::vector<std::thread> threads;

        // Create threads for processing the negative manifold
        for (int i = 0; i < num_threads; ++i) {
            int start_index = i * chunk_size_neg;
            int end_index = (i + 1) * chunk_size_neg;
            if (i == num_threads - 1) end_index = ICs_Stable_Negative.size(); // Handle remainder
            threads.push_back(std::thread(process_negative_manifold, start_index, end_index, std::ref(ICs_Stable_Negative),
                                          mu, poincare_section_x, std::ref(stable_negative_state), std::ref(hit_poincare_section)));
        }

        // Create threads for processing the positive manifold
        for (int i = 0; i < num_threads; ++i) {
            int start_index = i * chunk_size_pos;
            int end_index = (i + 1) * chunk_size_pos;
            if (i == num_threads - 1) end_index = ICs_Stable_Positive.size(); // Handle remainder
            threads.push_back(std::thread(process_positive_manifold, start_index, end_index, std::ref(ICs_Stable_Positive),
                                          mu, poincare_section_x, std::ref(stable_positive_state), std::ref(hit_poincare_section)));
        }

        // Join all threads
        for (auto& t : threads) {
            t.join();
        }
    
    
    std::vector<std::vector<double>> stable_manifold;
*/
    
    double stable_negative_state[number_of_manifolds][7] = {100};
    double stable_positive_state[number_of_manifolds][7] = {100};
    
        
    std::vector<std::vector<double>> stable_manifold;
    
    double poincare_section_x = -mu;
    int hit_poincare_section = 0;
    auto arrays_neg = Generate_Negative_Manifolds( Initial_Conditions,  Tp,  mu,  number_of_manifolds);
    std::cout << "The Initial Conditions for the negative portion of the manifold have been generated" << std::endl;
    std::vector<std::vector<double>> ICs_Stable_Negative = arrays_neg.first;
    std::vector<std::vector<double>> ICs_Unstable_Negative = arrays_neg.second;
    
    for (int number_of_negative_IC = 0; number_of_negative_IC < ICs_Stable_Negative.size(); ++number_of_negative_IC){
        ICs[0] = ICs_Stable_Negative[number_of_negative_IC][0];
        ICs[1] = ICs_Stable_Negative[number_of_negative_IC][1];
        ICs[2] = ICs_Stable_Negative[number_of_negative_IC][2];
        ICs[3] = ICs_Stable_Negative[number_of_negative_IC][3];
        ICs[4] = ICs_Stable_Negative[number_of_negative_IC][4];
        ICs[5] = ICs_Stable_Negative[number_of_negative_IC][5];
        std::vector<double> Output = Get_Poincare_Section_2( mu,ICs,poincare_section_x);

        if(Output[0] != 0){
            
            stable_negative_state[hit_poincare_section][0] = Output[0];
            stable_negative_state[hit_poincare_section][1] = Output[1];
            stable_negative_state[hit_poincare_section][2] = Output[2];
            stable_negative_state[hit_poincare_section][3] = Output[3];
            stable_negative_state[hit_poincare_section][4] = Output[4];
            stable_negative_state[hit_poincare_section][5] = Output[5];
            stable_negative_state[hit_poincare_section][6] = Output[6];
            hit_poincare_section = hit_poincare_section+1;
            
          //  stable_negative_state.push_back({Output[0],Output[1],Output[2],Output[3],Output[4],Output[5],Output[6],});
            
            
            Output.clear();
            
        }
        std::fill(ICs.begin(),ICs.end(),0);
        double percent = number_of_negative_IC/150.*100.0;
        cout << "The negative portion of the manifold is: " << percent <<"% finished" << endl;
    }
    
    
    hit_poincare_section = 0;
    auto arrays_pos = Generate_Positive_Manifolds( Initial_Conditions,  Tp,  mu,  number_of_manifolds);
    std::cout << "The Initial Conditions for the positive portion of the manifold have been generated" << std::endl;

    std::vector<std::vector<double>> ICs_Stable_Positive = arrays_pos.first;
    std::vector<std::vector<double>> ICs_Unstable_Positive = arrays_pos.second;
    
    for (int number_of_positive_IC = 0; number_of_positive_IC < ICs_Stable_Positive.size(); ++number_of_positive_IC){
        ICs[0] = ICs_Stable_Positive[number_of_positive_IC][0];
        ICs[1] = ICs_Stable_Positive[number_of_positive_IC][1];
        ICs[2] = ICs_Stable_Positive[number_of_positive_IC][2];
        ICs[3] = ICs_Stable_Positive[number_of_positive_IC][3];
        ICs[4] = ICs_Stable_Positive[number_of_positive_IC][4];
        ICs[5] = ICs_Stable_Positive[number_of_positive_IC][5];
        std::vector<double> Output = Get_Poincare_Section_2( mu,ICs,poincare_section_x);
        if(Output[0] != 0){
            stable_positive_state[hit_poincare_section][0] = Output[0];
            stable_positive_state[hit_poincare_section][1] = Output[1];
            stable_positive_state[hit_poincare_section][2] = Output[2];
            stable_positive_state[hit_poincare_section][3] = Output[3];
            stable_positive_state[hit_poincare_section][4] = Output[4];
            stable_positive_state[hit_poincare_section][5] = Output[5];
            stable_positive_state[hit_poincare_section][6] = Output[6];
            hit_poincare_section = hit_poincare_section+1;
            Output.clear();
        }
        std::fill(ICs.begin(),ICs.end(),0);
        double percent = number_of_positive_IC/150.*100.0;
        cout << "The positive portion of the manifold is: " << percent <<"% finished" << endl;
    }
     
     
    //]
    
    
    //[ Finding which manifold (positive or negative) ends closer to the Earth
    double y_pos_min = 101;
    double y_neg_min = 101;
    for( int i = 0; i < number_of_manifolds; ++i){
        if(stable_positive_state[i][1] != 0){
        if(std::abs(stable_positive_state[i][1]) < y_pos_min){
            y_pos_min = std::abs(stable_positive_state[i][1]);
        }
    }
    }
    for( int i = 0; i < number_of_manifolds; ++i){
        if(stable_negative_state[i][1] != 0){
            
            if(std::abs(stable_negative_state[i][1]) < y_neg_min){
                y_neg_min = std::abs(stable_negative_state[i][1]);
            }
        }
    }
   
    if( y_pos_min > y_neg_min){
        for(int number_of_stable_manifolds = 0; number_of_stable_manifolds < number_of_manifolds; number_of_stable_manifolds++){
            if(stable_negative_state[number_of_stable_manifolds][0] != 100){
                
                stable_manifold.push_back({
                    stable_negative_state[number_of_stable_manifolds][0],
                    stable_negative_state[number_of_stable_manifolds][1],
                    stable_negative_state[number_of_stable_manifolds][2],
                    stable_negative_state[number_of_stable_manifolds][3],
                    stable_negative_state[number_of_stable_manifolds][4],
                    stable_negative_state[number_of_stable_manifolds][5],
                    stable_negative_state[number_of_stable_manifolds][6]
                });
                T2.push_back(stable_positive_state[number_of_stable_manifolds][6]*hour);
                number_of_manifolds_at_poincare = number_of_manifolds_at_poincare+1;
            }
        }
        
    }else{
        for(int number_of_stable_manifolds = 0; number_of_stable_manifolds < number_of_manifolds; number_of_stable_manifolds++){
            if(stable_negative_state[number_of_stable_manifolds][0] != 100){
                
                stable_manifold.push_back({
                    stable_positive_state[number_of_stable_manifolds][0],
                    stable_positive_state[number_of_stable_manifolds][1],
                    stable_positive_state[number_of_stable_manifolds][2],
                    stable_positive_state[number_of_stable_manifolds][3],
                    stable_positive_state[number_of_stable_manifolds][4],
                    stable_positive_state[number_of_stable_manifolds][5],
                    stable_positive_state[number_of_stable_manifolds][6]
                });
                T2.push_back(stable_positive_state[number_of_stable_manifolds][6]*hour);
                number_of_manifolds_at_poincare = number_of_manifolds_at_poincare+1;
                
            }
        }
    }
    
    
    
    //]
    //[ Finding if any branches of the manifold intersect with the moon
/*
        int number_of_hits = 0;
    size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = stable_manifold.size() / num_threads;

        std::vector<std::thread> threads;

        // Launch threads
        for (size_t t = 0; t < num_threads; ++t) {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? stable_manifold.size() : start + chunk_size;

            threads.emplace_back(Check_If_Hits_Moon, std::cref(stable_manifold),
                                 std::ref(hits_moon), std::ref(number_of_hits),
                                 start, end, abs_err, rel_err, LU, mu);
        }

        // Join threads
        for (auto& thread : threads) {
            thread.join();
        }
    */
    
    
    double dt;
    double distance_moon;
    int number_of_hits = 0;
    for(int j = 0; j < stable_manifold.size(); ++j){
        std::fill(w.begin(),w.end(),0);
        w[0] = stable_manifold[j][0];
        w[1] = stable_manifold[j][1];
        w[2] = stable_manifold[j][2];
        w[3] = stable_manifold[j][3];
        w[4] = stable_manifold[j][4];
        w[5] = stable_manifold[j][5];
        TOF = stable_manifold[j][6];
        
        dt = TOF/4999.0;
        std::vector<double> times_vector(5000);
        for (size_t p = 0; p < 5000; ++p) {
            times_vector[p] = p * dt;
        }
        w_vec.clear();
        times.clear();
        
        integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop, w, times_vector.begin(),times_vector.end(), dt, push_back_state_and_time(w_vec, times));
        
     
        for( int moon_test = 0; moon_test <w_vec.size(); ++moon_test){
            distance_moon = pow(pow((w_vec[moon_test][0]-(1.0-mu)),2.0)+pow(w_vec[moon_test][1],2.0)+pow(w_vec[moon_test][2],2.0),0.5);
            
            if(distance_moon < 1740/LU){
                hits_moon.push_back(j);
                number_of_hits = number_of_hits+1;
                break;
            }
            
        }
    }
     
    //]
    
    
//[ Finding the initial guess for delta_V and TOF;
    
    std::vector<double>r_ECI = state1.r2_vec;
    std::vector<double>v_ECI = state1.v2_vec;
    
    std::vector<double> r_rot = ECI2Rotating_Frame_Rotation(r_ECI[0],r_ECI[1],r_ECI[2],v_ECI[0],v_ECI[1],v_ECI[2],0);
    std::vector<double> r_target = {r_rot[0], r_rot[1], r_rot[2]};
    std::vector<double> v_target = {r_rot[3], r_rot[4], r_rot[5]};

    
    
int first_manifold_branch = 1;
int good_first_guess = 1;



w[0] = Initial_Conditions[0];
w[1] = Initial_Conditions[1];
w[2] = Initial_Conditions[2];
w[3] = Initial_Conditions[3];
w[4] = Initial_Conditions[4];
w[5] = Initial_Conditions[5];
double dt_initial_guess = Tp/999.0;
std::vector<double> times_vector_initial_guess(1000);
for (size_t p = 0; p < 1000; ++p) {
    times_vector_initial_guess[p] = p * dt_initial_guess;
}
w_vec.clear();
times.clear();

integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop, w, times_vector_initial_guess.begin(),times_vector_initial_guess.end(), dt_initial_guess, push_back_state_and_time(w_vec, times));
int departure_point = 0;

double departure_x = w_vec[departure_point][0];
double departure_y = w_vec[departure_point][1];
double departure_z = w_vec[departure_point][2];

    
    

while(good_first_guess == 1){
    std::vector<double> IC0(6);
    IC0[0] = stable_manifold[first_manifold_branch][0];
    IC0[1] = stable_manifold[first_manifold_branch][1];
    IC0[2] = stable_manifold[first_manifold_branch][2];
    IC0[3] = stable_manifold[first_manifold_branch][3];
    IC0[4] = stable_manifold[first_manifold_branch][4];
    IC0[5] = stable_manifold[first_manifold_branch][5];
    
    double v0_man = pow(pow(IC0[3],2.0)+pow(IC0[4],2.0)+pow(IC0[5],2.0),0.5);
    double v0_target = pow(pow(v_target[0],2.0)+pow(v_target[1],2.0)+pow(v_target[2],2.0),0.5);

    double delta_V0 = abs(v0_target-v0_man);
    
    double TOF_low = 0.1;
    double TOF_high = 3.0;
    double delta_V_low = -(delta_V0 + delta_V0*2);
    double delta_V_high = delta_V0 + delta_V0*6;
    
    
    int number_of_points = 40;
    
    
    std::vector<double> Output = Get_DeltaV_TOF_Initial_Guess(r_rot, IC0,  delta_V_high,  delta_V_low,  TOF_high,  TOF_low, number_of_points, LU, mu);

    Tps = Output[0];
    TOF = Tps;
    v2_0 = Output[1];
    cout << Tps << endl;
    cout << v2_0 << endl;
    

    
    if(Output[0] == 0){
        first_manifold_branch = floor(number_of_manifolds/4)+first_manifold_branch;
        
    }
    else{
        good_first_guess = 0;
    }
    //]

    cout << "branch number: "<< first_manifold_branch << endl;
    
}
    
    
    
    
    //[ Differential Corrections
    
    std::vector<double> v(3);
    std::vector<double> r(3);
    std::vector<double> IC_transfer(6);
    std::vector<double> IC_X(6);
    v[0] = stable_manifold[0][3];
    v[1] = stable_manifold[0][4];
    v[2] = stable_manifold[0][5];
    double norm_v = std::sqrt(pow(v[0],2.0)+pow(v[1],2.0)+pow(v[2],2.0));
    double xdot0 = v2_0*v[0]/norm_v+v[0];
    double ydot0 = v2_0*v[1]/norm_v+v[1];
    double zdot0 = v2_0*v[2]/norm_v+v[2];
    double err = 1;
    int k = 1;
    int max_iterations = 40;
    bool PrintMetaData = 1;

for(int i = 0; i < stable_manifold.size(); ++i){
    std::fill(w.begin(),w.end(),0);
    w_vec.clear();
    times.clear();
    
    r[0] = stable_manifold[i][0];
    r[1] = stable_manifold[i][1];
    r[2] = stable_manifold[i][2];
    
    v[0] = stable_manifold[i][3];
    v[1] = stable_manifold[i][4];
    v[2] = stable_manifold[i][5];
    norm_v = std::sqrt(pow(v[0],2.0)+pow(v[1],2.0)+pow(v[2],2.0));
    
    if (std::sqrt(pow(r[0]-(1-mu),2.0)+pow(r[1],2.0)+pow(r[2],2.0)) < 1740.0/LU){
        cout << "started too close" << endl;
    continue;
}
        
   
    
    
    
        
        if (j == max_iterations || TOF <= 0 || err > 10){
            //cout << "Skipped this branch.  TOF:  " << TOF << ".  Error: " << err <<endl;
            xdot0 = v2_0*v[0]/norm_v+v[0];
            ydot0 = v2_0*v[1]/norm_v+v[1];
            zdot0 = v2_0*v[2]/norm_v+v[2];
            TOF = Tps;
            IC_X[0] = r[0];
            IC_X[1] = r[1];
            IC_X[2] = r[2];
            IC_X[3] = xdot0;
            IC_X[4] = ydot0;
            IC_X[5] = zdot0;
            err = 1;
        }
        IC_X[0] = r[0];
        IC_X[1] = r[1];
        IC_X[2] = r[2];
        IC_X[3] = xdot0;
        IC_X[4] = ydot0;
        IC_X[5] = zdot0;
        
        j = 1;
        while (j < max_iterations && TOF > 0 && err < 10){
            
            w[0] = IC_X[0];
            w[1] = IC_X[1];
            w[2] = IC_X[2];
            w[3] = IC_X[3];
            w[4] = IC_X[4];
            w[5] = IC_X[5];
            
            w[6] = w[13] = w[20] = w[27] = w[34] = w[41] = 1;
            w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[14] = w[15] = w[16] = w[17] = w[18] = w[19] = w[21] = w[22] = w[23] = w[24] = w[25] = w[26] = w[28] = w[29] = w[30] = w[31] = w[32] = w[33] = w[35] = w[36] = w[37] = w[38] = w[39] = w[40] = 0;
            
            
           
           
           
            //[ define_adapt_stepper
            typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
            auto controlled_stepper =make_controlled(abs_err,rel_err, error_stepper_type());
            //]
     
            
            
          
            double t_start = 0.0;
            size_t max_steps = 1000;  // Set the maximum number of steps
            double dt2 = TOF/(max_steps-1);
            
            std::vector<double> times_vector2(1000);
            for (size_t p = 0; p < 1000; ++p) {
                times_vector2[p] = p * dt2;
            }
            
            integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_STM_prop_backwards, w, times_vector2.begin(),times_vector2.end(), dt2, push_back_state_and_time(w_vec, times));

            
            //[ Finding the State Transition Matrix and State at the end of the integration
            int end_of_int = w_vec.size();
            end_of_int = end_of_int -1;
                
            ublas::matrix<double>PHI_final(6,6);
            PHI_final(0,0) = w_vec[end_of_int][6]; PHI_final(0,1) = w_vec[end_of_int][7]; PHI_final(0,2) = w_vec[end_of_int][8]; PHI_final(0,3) = w_vec[end_of_int][9]; PHI_final(0,4) = w_vec[end_of_int][10]; PHI_final(0,5) = w_vec[end_of_int][11];
            PHI_final(1,0) = w_vec[end_of_int][12]; PHI_final(1,1) = w_vec[end_of_int][13]; PHI_final(1,2) = w_vec[end_of_int][14]; PHI_final(1,3) = w_vec[end_of_int][15]; PHI_final(1,4) = w_vec[end_of_int][16]; PHI_final(1,5) = w_vec[end_of_int][17];
            PHI_final(2,0) = w_vec[end_of_int][18]; PHI_final(2,1) = w_vec[end_of_int][19]; PHI_final(2,2) = w_vec[end_of_int][20]; PHI_final(2,3) = w_vec[end_of_int][21]; PHI_final(2,4) = w_vec[end_of_int][22]; PHI_final(2,5) = w_vec[end_of_int][23];
            PHI_final(3,0) = w_vec[end_of_int][24]; PHI_final(3,1) = w_vec[end_of_int][25]; PHI_final(3,2) = w_vec[end_of_int][26]; PHI_final(3,3) = w_vec[end_of_int][27]; PHI_final(3,4) = w_vec[end_of_int][28]; PHI_final(3,5) = w_vec[end_of_int][29];
            PHI_final(4,0) = w_vec[end_of_int][30]; PHI_final(4,1) = w_vec[end_of_int][31]; PHI_final(4,2) = w_vec[end_of_int][32]; PHI_final(4,3) = w_vec[end_of_int][33]; PHI_final(4,4) = w_vec[end_of_int][34]; PHI_final(4,5) = w_vec[end_of_int][35];
            PHI_final(5,0) = w_vec[end_of_int][36]; PHI_final(5,1) = w_vec[end_of_int][37]; PHI_final(5,2) = w_vec[end_of_int][38]; PHI_final(5,3) = w_vec[end_of_int][39]; PHI_final(5,4) = w_vec[end_of_int][40]; PHI_final(5,5) = w_vec[end_of_int][41];
            double x = w_vec[end_of_int][0];
            double y = w_vec[end_of_int][1];
            double z = w_vec[end_of_int][2];
            double xdot = w_vec[end_of_int][3];
            double ydot = w_vec[end_of_int][4];
            double zdot = w_vec[end_of_int][5];
                
                
                
            std::vector<double> state0 (6);
            state0[0] = w_vec[0][0];
            state0[1] = w_vec[0][1];
            state0[2] = w_vec[0][2];
            state0[3] = w_vec[0][3];
            state0[4] = w_vec[0][4];
            state0[5] = w_vec[0][5];
            double v2 = std::sqrt(pow(state0[3],2.0)+pow(state0[4],2.0)+pow(state0[5],2.0));
            w_vec.clear();
            times.clear();
            //]
            
            double p1 = pow(pow(x+mu,2)+pow(y,2)+pow(z,2),0.5);
            double p2 = pow(pow(x+mu-1,2)+pow(y,2)+pow(z,2),0.5);
            
            vector< double > vector_field(6);
            vector_field[0] = xdot;
            vector_field[1] = ydot;
            vector_field[2] = zdot;
            vector_field[3] = 2*ydot+x -(1-mu)*(x+mu)/pow(p1,3)-(mu)*(x-1+mu)/pow(p2,3);
            vector_field[4] = -2*xdot+y-(1-mu)*(y)/pow(p1,3)-(mu)*(y)/pow(p2,3);
            vector_field[5] = -(1-mu)*(z)/pow(p1,3)-(mu)*(z)/pow(p2,3);
            
            
            Eigen::MatrixXd Augmented_STM(3,4);
            Augmented_STM(0,0) = PHI_final(0,3);
            Augmented_STM(0,1) = PHI_final(0,4);
            Augmented_STM(0,1) = PHI_final(0,5);
            Augmented_STM(0,3) = xdot;
            
            Augmented_STM(1,0) = PHI_final(1,3);
            Augmented_STM(1,1) = PHI_final(1,4);
            Augmented_STM(1,2) = PHI_final(1,5);
            Augmented_STM(1,3) = ydot;
             
            Augmented_STM(2,0) = PHI_final(2,3);
            Augmented_STM(2,1) = PHI_final(2,4);
            Augmented_STM(2,2) = PHI_final(2,5);
            Augmented_STM(2,3) = zdot;
            
            Eigen::MatrixXd pinv = pseudoInverse(Augmented_STM);
            ublas::matrix<double>Constraint(3,1);
            Constraint(0,0) = x-r_target[0];
            Constraint(1,0) = y-r_target[1];
            Constraint(2,0) = z-r_target[2];


            err = std::sqrt(pow(Constraint(0,0),2)+pow(Constraint(1,0),2)+pow(Constraint(2,0),2));
            ublas::matrix<double> inverse_STM(4,3);
            
            
            inverse_STM(0,0) = pinv(0,0);
            inverse_STM(0,1) = pinv(0,1);
            inverse_STM(0,2) = pinv(0,2);
            
            inverse_STM(1,0) = pinv(1,0);
            inverse_STM(1,1) = pinv(1,1);
            inverse_STM(1,2) = pinv(1,2);

            inverse_STM(2,0) = pinv(2,0);
            inverse_STM(2,1) = pinv(2,1);
            inverse_STM(2,2) = pinv(2,2);

            inverse_STM(3,0) = pinv(3,0);
            inverse_STM(3,1) = pinv(3,1);
            inverse_STM(3,2) = pinv(3,2);
            
            ublas::matrix<double>DF(4,1);
            DF = ublas::prod(inverse_STM,Constraint);
            

            if (std::abs(err) < tol) {
                w_vec2.clear();
                std::fill(w.begin(),w.end(),0);
                w[0] = IC_X[0];
                w[1] = IC_X[1];
                w[2] = IC_X[2];
                w[3] = IC_X[3];
                w[4] = IC_X[4];
                w[5] = IC_X[5];
                
                w[6] = w[13] = w[20] = w[27] = w[34] = w[41] = 1;
                w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[14] = w[15] = w[16] = w[17] = w[18] = w[19] = w[21] = w[22] = w[23] = w[24] = w[25] = w[26] = w[28] = w[29] = w[30] = w[31] = w[32] = w[33] = w[35] = w[36] = w[37] = w[38] = w[39] = w[40] = 0;
                double dt3 = TOF/2999;
                std::vector<double> times_vector3(3000);
                for (size_t p = 0; p < 3000; ++p) {
                    times_vector3[p] = p * dt3;
                }
                //[ Integration to find the STM at one full time period
                integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_prop_backwards, w,  times_vector3.begin(),times_vector3.end(), dt3, push_back_state_and_time(w_vec3, times));
                //]
                
                int hits_earth_moon = 0;
                double distance_Earth = 0;
                for(int test = 0; test < w_vec3.size(); ++test){
                    distance_Earth = std::sqrt(pow((w_vec3[test][0]-(-mu)),2.0)+pow(w_vec3[test][1],2.0)+pow(w_vec3[test][2],2.0));
                    if(distance_Earth < 6400/LU || std::count(hits_moon.begin(),hits_moon.end(), i) > 0) {
                        delta_vx[k] = 0;
                        delta_vy[k] = 0;
                        delta_vz[k] = 0;
                        
                        x_departure[k] = 0;
                        y_departure[k] = 0;
                        z_departure[k] = 0;
                        xdot_departure[k] = 0;
                        ydot_departure[k] = 0;
                        zdot_departure[k] = 0;
                        T1[k] = 0;
                        T2_elapsed[k] = 0;
                        hits_earth_moon = 1;
                     //   cout << "Hits the Earth" << endl;
                        break;
                    }
                    
                }
                
                
                if(hits_earth_moon == 0 && TOF > 0 ){
                    delta_vx.push_back(v[0]-IC_X[3]);
                    delta_vy.push_back(v[1]-IC_X[4]);
                    delta_vz.push_back(v[2]-IC_X[5]);
                    
                    x_departure.push_back(x);
                    y_departure.push_back(y);
                    z_departure.push_back(z);
                    xdot_departure.push_back(xdot);
                    ydot_departure.push_back(ydot);
                    zdot_departure.push_back(zdot);
                    T1.push_back(TOF*hour);
                    T2_elapsed.push_back(T1[k] + T2[i]);
                    
                    
                    if(std::sqrt(pow(delta_vx[k],2)+pow(delta_vy[k],2)+pow(delta_vz[k],2)) < 5 ){
                        //cout << std::fixed << std::setprecision(14) << x_departure.back() << '\t' << y_departure.back() << '\t' <<  z_departure.back() << '\t' << xdot_departure.back() << '\t' << ydot_departure.back() << '\t' << zdot_departure.back() << '\t' << T1.back() << '\t' << delta_vx.back() << '\t' << delta_vy.back() << '\t' << delta_vz.back() << '\t' << T2_elapsed.back() << '\n';
                        
                        //outputFile << std::fixed << std::setprecision(15) << x_departure[k] << '\t' << y_departure[k] << '\t' <<  z_departure[k] << '\t' << xdot_departure[k] << '\t' << ydot_departure[k] << '\t' << zdot_departure[k] << '\t' << T1[k] << '\t' << delta_vx[k] << '\t' << delta_vy[k] << '\t' << delta_vz[k] << '\t' << T2_elapsed[k] << '\n';
                        Create_Outputs( state1, arrival_orbit,  LU,  TU,  mu, x_departure.back(), y_departure.back(), z_departure.back(), xdot_departure.back(), ydot_departure.back(), zdot_departure.back(), T1.back(), delta_vx.back(), delta_vy.back(), delta_vz.back(), T2_elapsed.back(), PrintMetaData);
                        PrintMetaData = 0;
                    }
                    //]
                }
                
                
                k = k+1;
                
                
                break;
                
            } // end of the outputs for loop
            
            

            //[ Initial condition and time period corrections
          
            IC_X[3] = IC_X[3] - 0.95*DF(0,0);
            IC_X[4] = IC_X[4] - 0.95*DF(1,0);
            IC_X[5] = IC_X[5] - 0.95*DF(2,0);
            TOF = TOF + 0.95*DF(3,0);
            //]
            
            cout << "Iteration number:" << j <<  "   err: " << err << endl;
            j = j+1;
            
       
            
        } //end of each correction loop
        
        cout << "Manifold Branch Number: " << i << endl;

      
        
        
    }//end of loop first delta_z family loop
    
    
    
    
    
    
    
    
    
    
} //end
