//
//  Proccessing_Negative_Manifold.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/12/24.
//

#ifndef Proccessing_Negative_Manifold_hpp
#define Proccessing_Negative_Manifold_hpp

#include <vector>
#include <thread>
#include <mutex>
#include <iostream>
#include <algorithm>
void process_negative_manifold(int start_index, int end_index, const std::vector<std::vector<double>>& ICs_Stable_Negative,
                               double mu, double poincare_section_x, std::vector<std::vector<double>>& stable_negative_state);
#endif /* Proccessing_Negative_Manifold_hpp */
