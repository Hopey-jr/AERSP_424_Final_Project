//
//  Transfers.hpp
//  Final_Project_V2
//
//  Created by jonathon hope on 12/10/24.
//

#ifndef Transfers_hpp
#define Transfers_hpp

#include "IOD.hpp"          // Include the IOD class header
#include "Arrival_Orbit.hpp" // Include the Arrival_Orbit class header
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



// Function declaration
void Create_Transfers(IOD& state1, Arrival_Orbit& arrival_orbit);
#endif /* Transfers_hpp */
