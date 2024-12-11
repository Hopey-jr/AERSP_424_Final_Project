//
//  EOM_function.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#include "EOM_function.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>


typedef std::vector< double > state_type;
namespace ublas = boost::numeric::ublas;

void state_prop( const state_type &w , state_type &dwdt , const double /* t */ )
{
     double mu = 1.215058560962404E-2;

    dwdt[0] = w[3];
    dwdt[1] = w[4];
    dwdt[2] = w[5];
 
    double P1;
    double P2;
    double Uy;
    double Ux;
     P1=pow(pow(w[0]+mu,2)+pow(w[1],2)+pow(w[2],2),.5);
     P2=pow(pow(w[0]-1+mu,2)+pow(w[1],2)+pow(w[2],2),.5);
     Ux=w[0]-((1-mu)*(w[0]+mu))/pow(P1,3)-(mu*(w[0]-1+mu))/pow(P2,3);
     Uy=w[1]-((1-mu)*w[1])/pow(P1,3)-mu*w[1]/pow(P2,3);
    dwdt[3] =  2*w[4]+Ux;
    dwdt[4] =  -2*w[3]+Uy;
    dwdt[5] = -(1-mu)*w[2]/pow(P1,3)-mu*w[2]/pow(P2,3);
    
    
    
}
