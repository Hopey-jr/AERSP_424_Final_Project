//
//  Hits_Moon_Check.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/11/24.
//

#ifndef Hits_Moon_Check_hpp
#define Hits_Moon_Check_hpp

#include <stdio.h>
#include <vector>
#include <mutex>
#include <thread>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "PushBackStateAndTime.hpp"
#include "EOM_function.hpp"

void Check_If_Hits_Moon( std::vector<std::vector<double>>& stable_manifold, std::vector<double>& hits_moon, int start, int end, double abs_err, double rel_err, double LU, double mu);
#endif /* Hits_Moon_Check_hpp */
