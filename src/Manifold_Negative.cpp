//  Description:
//               The purpose of this code is to determine the initial
//               conditions for both the stable and unstable branches of the
//               negative manifold. The stable part of the manifold will go
//               towards the periodic orbit and the unstable part of the
//               manifold will go away from the periodic orbit


//Inputs:
//       Initial_Conditions - The initial conditions for the periodic orbit coming from the
//       orbit catalog
//       Tp - The orbital period of the periodic orbit coming from the
//       orbit catalog
//       mu - The nondimensional reduced mass of the Earth-Moon system [LU^3/TU^2]
//       n - The number of points the user chooses to discretize the
//       periodic orbit. This number will also be the number of branches in
//       the manifold
//Output:
//       ICs_stable - The initial conditions for the stable portion of the
//       manifold
//       ICs_unstable - The initial conditions for the unstable portion of the
//       manifold
//Locals:
//       LU - The characteristic length unit derived from the average distance between the Earth and Moon
//       w_vec - The state and value of the STM from 0 to tau
//       w - The vector and variable used to hold the initital conditions for the integration that include both the initital state and initial STM
//       times - The timestamps of each integration
//       dt - The initial timestep used for the integration from 0 to Tp
//       times_vector - The time vector used to propagate the periodic orbit
//       end_of_int - Used to denote the final state and STM after the integration
//       M - The monodromy matrix
//       eigenvectors - The eigenvectors of the monodromy matrix
//       eigenvalues - The eigenvalues of the monodromy matrix
//       min_val - The minimum value of the eigenvalues
//       max_val - The maximum value of the eigenvalues
//       inner - The position of the eigenvalue that corresponds to the stable eigenvalue and manifold
//       outer - The position of the eigenvalue that corresponds to the unstable eigenvalue and manifold
//       eigenvector_stable - The eigenvector corresponding to the stable manifold
//       eigenvector_unstable - The eigenvector corresponding to the unstable manifold
//       x_vec - The vector that holds all the x-positions from the propagation of the periodic orbit
//       y_vec - The vector that holds all the y-positions from the propagation of the periodic orbit
//       z_vec - The vector that holds all the z-positions from the propagation of the periodic orbit
//       xdot_vec - The vector that holds all the x-component of velocity from the propagation of the periodic orbit
//       ydot_vec - The vector that holds all the y-component of velocity from the propagation of the periodic orbit
//       zdot_vec - The vector that holds all the z-component of velocity from the propagation of the periodic orbit
//       time_vec - The vector that holds all the timestamps from the propagation of the periodic orbit
//       epsilon - The scaling factor of the perturbation used to create the initial conditions
//       IC_x - The x-position at a given timestamp around the periodic orbit
//       IC_y - The y-position at a given timestamp around the periodic orbit
//       IC_z - The z-position at a given timestamp around the periodic orbit
//       IC_xdot - The x-component of velocity at a given timestamp around the periodic orbit
//       IC_ydot - The y-component of velocity at a given timestamp around the periodic orbit
//       IC_zdot - The z-component of velocity at a given timestamp around the periodic orbit
//       T_PO - The time at a given timestamp around the periodic orbit
//       dt2 - The timestep from 0 to T_PO along the periodic orbit used to find the initial conditions of each branch of the manifold
//       times_vector2 - The time vector used to find the initial conditions of each branch of the manifold
//       STM - The state transition matrix from 0 to 'end_time'
//       pert_stable - The  perturbation used to create the initial conditions for the stable manifold
//       pert_unstable - The  perturbation used to create the initial conditions for the unstable manifold
//       norm_pert_stable - The magnitude of the stable perturbation
//       norm_pert_unstable - The magnitude of the unstable perturbation

//Dependent Functions:
//       EOM_STM_function.h" - Propagates the STM and state from t0 to tf

//function call:
// auto arrays = Generate_Negative_Manifolds( Initial_Conditions,  Tp,  mu,  n);

#include "Manifold_Negative.hpp"
#include <cmath>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Eigenvalues>

#include "EOM_STM_function.h"
#include "EOM_STM_backwards_function.hpp"
#include "PushBackStateAndTime.hpp"



namespace ublas = boost::numeric::ublas;
using namespace boost::numeric::odeint;

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
Generate_Negative_Manifolds(const std::vector<double>& x0, double Tp, double mu, int n) {
    using namespace std;
   
    using namespace boost::numeric::odeint;
    double abs_err = 1.0e-13; //absolute tolerance
    double rel_err = 1.0e-13; //relative tolerance
    double LU = 384400;
    
    
    vector<state_type> w_vec;
    vector<state_type> w_vec2;
    
    vector<double> times;
    //[ state_initialization
    std::vector<std::vector<double>> ICs_stable(n, std::vector<double>(7,0.0));
    std::vector<std::vector<double>> ICs_unstable(n, std::vector<double>(7,0.0));

    
    double x_IC = x0[0];
    double y_IC = x0[1];
    double z_IC = x0[2];
    double xdot_IC = x0[3];
    double ydot_IC = x0[4];
    double zdot_IC = x0[5];
    
    state_type w(42);
    std::fill(w.begin(),w.end(),0);
    
    w[0] = x_IC;
    w[1] = y_IC;
    w[2] = z_IC;
    w[3] = xdot_IC;
    w[4] = ydot_IC;
    w[5] = zdot_IC;


    
    w[6] = w[13] = w[20] = w[27] = w[34] = w[41] = 1;
    w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[14] = w[15] = w[16] = w[17] = w[18] = w[19] = w[21] = w[22] = w[23] = w[24] = w[25] = w[26] = w[28] = w[29] = w[30] = w[31] = w[32] = w[33] = w[35] = w[36] = w[37] = w[38] = w[39] = w[40] = 0;
    
    
    double dt = Tp/(n-1);
    
    
    
    //[ define_adapt_stepper
    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    //]
    
    //[
    
    std::vector<double> times_vector(n);
    for (size_t m = 0; m < times_vector.size(); ++m) {
        times_vector[m] = m * dt;
    }
    
    integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_STM_prop, w, times_vector.begin(),times_vector.end(), dt, push_back_state_and_time(w_vec, times));
    
    
    //]
    
    //[ Finding the Monodromy Matrix and stable and unstable eigenvectors
    int end_of_int = w_vec.size();
    end_of_int = end_of_int -1;
   
    
    Eigen::MatrixXd M(6, 6);
    M <<  w_vec[end_of_int][6],  w_vec[end_of_int][7],  w_vec[end_of_int][8],  w_vec[end_of_int][9],  w_vec[end_of_int][10],  w_vec[end_of_int][11],
    w_vec[end_of_int][12], w_vec[end_of_int][13], w_vec[end_of_int][14], w_vec[end_of_int][15], w_vec[end_of_int][16], w_vec[end_of_int][17],
    w_vec[end_of_int][18], w_vec[end_of_int][19], w_vec[end_of_int][20], w_vec[end_of_int][21], w_vec[end_of_int][22], w_vec[end_of_int][23],
    w_vec[end_of_int][24], w_vec[end_of_int][25], w_vec[end_of_int][26], w_vec[end_of_int][27], w_vec[end_of_int][28], w_vec[end_of_int][29],
    w_vec[end_of_int][30], w_vec[end_of_int][31], w_vec[end_of_int][32], w_vec[end_of_int][33], w_vec[end_of_int][34], w_vec[end_of_int][35],
    w_vec[end_of_int][36], w_vec[end_of_int][37], w_vec[end_of_int][38], w_vec[end_of_int][39], w_vec[end_of_int][40], w_vec[end_of_int][41];
    
    Eigen::EigenSolver<Eigen::MatrixXd> solver(M);
    Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
    
    
 
    
    

    int inner = 0;
    int outer = 0;
    
    double min_val = std::abs(eigenvalues(0));
    double max_val = std::abs(eigenvalues(0));
    
    for (int i = 1; i < 6; i ++){
        
        if ( std::abs(eigenvalues(i)) < min_val){
            min_val = eigenvalues(i);
            inner = i;
        }
        if ( std::abs(eigenvalues(i)) > max_val){
            max_val = eigenvalues(i);
            outer = i;
        }
        
    }
    

    Eigen::VectorXd eigenvector_stable(6);
    eigenvector_stable << -1*eigenvectors(0,inner), -1*eigenvectors(1,inner), -1*eigenvectors(2,inner), -1*eigenvectors(3,inner), -1*eigenvectors(4,inner), -1*eigenvectors(5,inner);

    Eigen::VectorXd eigenvector_unstable(6,1);
    eigenvector_unstable << eigenvectors(0,outer), eigenvectors(1,outer), eigenvectors(2,outer), eigenvectors(3,outer), eigenvectors(4,outer), eigenvectors(5,outer);
    
    //]

    
    //[ Creating the arrival orbit

    std::vector<double> x_vec(n);
    std::vector<double> y_vec(n);
    std::vector<double> z_vec(n);
    std::vector<double> xdot_vec(n);
    std::vector<double> ydot_vec(n);
    std::vector<double> zdot_vec(n);
    std::vector<double> time_vec(n);
    
    
    for(int i = 0; i < n; i++){
        x_vec[i] =w_vec[i][0];
        y_vec[i] =w_vec[i][1];
        z_vec[i] =w_vec[i][2];
        xdot_vec[i] =w_vec[i][3];
        ydot_vec[i] =w_vec[i][4];
        zdot_vec[i] =w_vec[i][5];
        time_vec[i] = times[i];
    }
    
    //]
    
    double epsilon = 70/LU;
    

    ICs_stable[0][0] = x_vec[0] - epsilon*eigenvector_stable(0);
    ICs_stable[0][1] = y_vec[0] - epsilon*eigenvector_stable(1);
    ICs_stable[0][2] = z_vec[0] - epsilon*eigenvector_stable(2);
    ICs_stable[0][3] = xdot_vec[0] - epsilon*eigenvector_stable(3);
    ICs_stable[0][4] = ydot_vec[0] - epsilon*eigenvector_stable(4);
    ICs_stable[0][5] = zdot_vec[0] - epsilon*eigenvector_stable(5);

    ICs_unstable[0][0] = x_vec[0] - epsilon*eigenvector_unstable(0);
    ICs_unstable[0][1] = y_vec[0] - epsilon*eigenvector_unstable(1);
    ICs_unstable[0][2] = z_vec[0] - epsilon*eigenvector_unstable(2);
    ICs_unstable[0][3] = xdot_vec[0] - epsilon*eigenvector_unstable(3);
    ICs_unstable[0][4] = ydot_vec[0] - epsilon*eigenvector_unstable(4);
    ICs_unstable[0][5] = zdot_vec[0] - epsilon*eigenvector_unstable(5);
    
    
    
    for(int i = 1; i < n; i++){ //start of the for loop for the IC of each branch of the manifold
        std::fill(w.begin(),w.end(),0);
        w_vec.clear();
        times.clear();
        double IC_x =  x_IC;
        double IC_y =  y_IC;
        double IC_z =  z_IC;
        double IC_xdot =  xdot_IC;
        double IC_ydot = ydot_IC;
        double IC_zdot =  zdot_IC;
        double T_PO = time_vec[i];
        w[0] = IC_x;
        w[1] = IC_y;
        w[2] = IC_z;
        w[3] = IC_xdot;
        w[4] = IC_ydot;
        w[5] = IC_zdot;
        
        
        w[6] = w[13] = w[20] = w[27] = w[34] = w[41] = 1;
        w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[14] = w[15] = w[16] = w[17] = w[18] = w[19] = w[21] = w[22] = w[23] = w[24] = w[25] = w[26] = w[28] = w[29] = w[30] = w[31] = w[32] = w[33] = w[35] = w[36] = w[37] = w[38] = w[39] = w[40] = 0;
        
        double dt2 = T_PO/999;
        std::vector<double> times_vector2(1000);
        for (size_t m = 0; m < times_vector2.size(); ++m) {
            times_vector2[m] = m * dt2;
        }

            integrate_times(make_controlled(abs_err, rel_err, error_stepper_type()), state_STM_prop, w, times_vector2.begin(),times_vector2.end(), dt2, push_back_state_and_time(w_vec, times));
            
            end_of_int = w_vec.size();
            end_of_int = end_of_int -1;


            Eigen::MatrixXd STM(6, 6);
            STM <<  w_vec[end_of_int][6],  w_vec[end_of_int][7],  w_vec[end_of_int][8],  w_vec[end_of_int][9],  w_vec[end_of_int][10],  w_vec[end_of_int][11],
            w_vec[end_of_int][12], w_vec[end_of_int][13], w_vec[end_of_int][14], w_vec[end_of_int][15], w_vec[end_of_int][16], w_vec[end_of_int][17],
            w_vec[end_of_int][18], w_vec[end_of_int][19], w_vec[end_of_int][20], w_vec[end_of_int][21], w_vec[end_of_int][22], w_vec[end_of_int][23],
            w_vec[end_of_int][24], w_vec[end_of_int][25], w_vec[end_of_int][26], w_vec[end_of_int][27], w_vec[end_of_int][28], w_vec[end_of_int][29],
            w_vec[end_of_int][30], w_vec[end_of_int][31], w_vec[end_of_int][32], w_vec[end_of_int][33], w_vec[end_of_int][34], w_vec[end_of_int][35],
            w_vec[end_of_int][36], w_vec[end_of_int][37], w_vec[end_of_int][38], w_vec[end_of_int][39], w_vec[end_of_int][40], w_vec[end_of_int][41];
            Eigen::VectorXd pert_stable = STM*eigenvector_stable;
        Eigen::VectorXd pert_unstable = STM*eigenvector_unstable;
        double norm_pert_stable = pert_stable.norm();
        pert_stable = pert_stable/norm_pert_stable;
        
        double norm_pert_unstable = pert_unstable.norm();
        pert_unstable = pert_unstable/norm_pert_unstable;
      

        
            ICs_stable[i][0] = x_vec[i] - epsilon*pert_stable(0);
            ICs_stable[i][1] = y_vec[i] - epsilon*pert_stable(1);
            ICs_stable[i][2] = z_vec[i] - epsilon*pert_stable(2);
            ICs_stable[i][3] = xdot_vec[i] - epsilon*pert_stable(3);
            ICs_stable[i][4] = ydot_vec[i] - epsilon*pert_stable(4);
            ICs_stable[i][5] = zdot_vec[i] - epsilon*pert_stable(5);
        
            ICs_unstable[i][0] = x_vec[i] - epsilon*pert_unstable(0);
            ICs_unstable[i][1] = y_vec[i] - epsilon*pert_unstable(1);
            ICs_unstable[i][2] = z_vec[i] - epsilon*pert_unstable(2);
            ICs_unstable[i][3] = xdot_vec[i] - epsilon*pert_unstable(3);
            ICs_unstable[i][4] = ydot_vec[i] - epsilon*pert_unstable(4);
            ICs_unstable[i][5] = zdot_vec[i] - epsilon*pert_unstable(5);
        
        }
        
    
    
    return std::make_pair(ICs_stable,ICs_unstable);
}
