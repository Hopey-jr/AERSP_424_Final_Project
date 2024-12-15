//
//  Orbit_Plot.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#include "Orbit_Plot.hpp"
#include <matplot/matplot.h>

#include <vector>
#include "Plot_Earth_Moon_Rotating.hpp"

void plot_position(std::vector<std::vector<double>>& x_vectors_all, std::vector<std::vector<double>>& y_vectors_all, std::vector<std::vector<double>>& z_vectors_all)
{
    using namespace matplot;
    hold("on");
    for(int i = 0; i < x_vectors_all.size(); ++i){
        plot3(x_vectors_all[i], y_vectors_all[i],z_vectors_all[i]);
     
    }
    Plot_earth_and_moon();
    // Formatting
//    title("Transfer from the spacecraft location to the arrival orbit");
    xlabel("x [LU]");
    ylabel("y [LU]");
    zlabel("z [LU]");
    hold("off");
    show();
    save("Transfer_Figure_X", "jpeg");
}


/*
 void plot_position(std::vector<double> x_vec, std::vector<double> y_vec, std::vector<double> z_vec)
 {
     using namespace matplot;
    

     plot3(x_vec, y_vec,z_vec);
     hold("on");
     

     // Formatting
 //    title("Transfer from the spacecraft location to the arrival orbit");
     xlabel("x [LU]");
     ylabel("y [LU]");
     zlabel("z [LU]");
     //hold("off");
     //show();
     //save("figure", "jpeg");
 }

 */
