//
//  ECI2Rotating_Frame_Coord_Transformation.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#ifndef ECI2Rotating_Frame_Coord_Transformation_hpp
#define ECI2Rotating_Frame_Coord_Transformation_hpp

#include <vector>
#include <cmath>

std::vector<double>ECI2Rotating_Frame_Rotation(double x_ECI, double y_ECI ,double z_ECI, double xdot_ECI, double ydot_ECI ,double zdot_ECI, double t);
#endif /* ECI2Rotating_Frame_Coord_Transformation_hpp */
