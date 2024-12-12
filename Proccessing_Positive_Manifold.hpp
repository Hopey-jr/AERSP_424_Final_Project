//
//  Proccessing_Positive_Manifold.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/12/24.
//

#ifndef Proccessing_Positive_Manifold_hpp
#define Proccessing_Positive_Manifold_hpp

#include <vector>
#include <thread>
#include <mutex>
#include <iostream>
#include <algorithm>
void process_positive_manifold(int start_index, int end_index, const std::vector<std::vector<double>>& ICs_Stable_Positive,
                               double mu, double poincare_section_x, std::vector<std::vector<double>>& stable_positive_state, int& hit_poincare_section);
#endif /* Proccessing_Positive_Manifold_hpp */
