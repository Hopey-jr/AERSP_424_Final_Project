//
//  Proccessing_Negative_Manifold.cpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/12/24.
//

#include "Proccessing_Negative_Manifold.hpp"
#include <vector>
#include <thread>
#include <mutex>
#include <iostream>
#include <algorithm>
#include "Poincare_Map_DC.hpp"

std::mutex mtx; // Mutex for thread-safe access to the stable states

void process_negative_manifold(int start_index, int end_index, const std::vector<std::vector<double>>& ICs_Stable_Negative,
                                double mu, double poincare_section_x, std::vector<std::vector<double>>& stable_negative_state, int& hit_poincare_section) {
    std::vector<double> ICs(6);
    
    for (int number_of_negative_IC = start_index; number_of_negative_IC < end_index; ++number_of_negative_IC) {
        ICs = ICs_Stable_Negative[number_of_negative_IC];
        std::vector<double> Output = Get_Poincare_Section_2(mu, ICs, poincare_section_x);

        if (Output[0] != 0) {
            std::lock_guard<std::mutex> lock(mtx); // Lock for thread-safe writing
            stable_negative_state[hit_poincare_section][0] = Output[0];
            stable_negative_state[hit_poincare_section][1] = Output[1];
            stable_negative_state[hit_poincare_section][2] = Output[2];
            stable_negative_state[hit_poincare_section][3] = Output[3];
            stable_negative_state[hit_poincare_section][4] = Output[4];
            stable_negative_state[hit_poincare_section][5] = Output[5];
            stable_negative_state[hit_poincare_section][6] = Output[6];
            hit_poincare_section++;
            Output.clear();
        }
    }
}
