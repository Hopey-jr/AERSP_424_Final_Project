# AERSP_424_Final_Project
NOTE FOR ANY GRADERS: Although the code will work for any arrival orbit in the catalog, there is no guarantee that any orbit will result in successfull transfers within the bounds of the problem. If any grader wants to test the full extent of the project, the following inputs will result in transfers that are catalogued and visualized: Measurement input = 1 | Orbit_Type = LYAPUNOV_L1 | Orbit_Number = 40. Warning: The runtime for these conditions is about 15 minutes. 


Description:
This purpose of this function is to design transfers from a spacecraft in a Earth-bound orbit to a Cislunar orbit around the L1 and L2 points. The function will begin by prompting the user for the following inputs: file name for the catalog output; measurement number based on text files prestored; arrival orbit type based on Lyapunov, Northern Halo, or Southern Halo at either the L1 or L2 point; and the orbit number based on how many orbits are in a prestored catalog. Once the user has all their inputs set, the code with start by finding the position and velocity based on the measurement set. After the program will build a manifold to ensure a ballistic transfer from the plane of the Earth all the way through the arrival orbit. The program will then use a differential correction scheme to find a transfer that goes from the end of the manifold to the location of the spacectaft. This ensures a consisent, low-energy transfer from the current location of the satellite to the chosen arrival orbit. Finally, each converged solution will be saved into a catalog and plotted using Matplot++ if the delta_v consumption is under 10 LU/TU. Throughout the code a logger will be used to provide checkpoints to ensure the code is running smoothly, as well as important outputs such as for the manifold and the delta_v consuption for a converged solution. The logger will also save key pieces of information to an external log. This project utilizes the Boost, Eigen, and Matplot++ libraries.


Dependent Functions: 

-Main Function:
The Main function will prompt the user to input the measurements and their choice of arrival orbit. This function will utilize classes to store the data. It will also call the function to generate the transfers

--Measurement_Input
This function assures that the user inputs a valid value for the measurement set

-- Orbit_Type_Input:
This function assures that the user inputs a valid value for the orbit type

-- Orbit_Number_Input:
This function assures that the user inputs a valid value for the orbit number

-IOD Function class:
This class will take the user input of the measurement set and create an object with the position and velocity of the spacecraft as the attributes.

--IOD_Measurement Function:
Collecting the terms from the user for the initial spacecraft measurements. For simplicity, the measurements are stored in text files and the user can specify which measurements   they want. The attributes of azimuth, elevation, position of the viewing site, and time of the viewing will be attached to an object within this class.
  
--Gauss Function:
This function will intake the measurements to solve for a distance to the spacecraft along with additional terms needed for the rest of the calculations. This function will also   need an additional Root_Finding function
  
--Root_Finding Function (Used online sources to help since the Boost and Eigen built-in functions didn’t work): 
This function utilizes Newton-Raphson iterations to solve for the real positive root of the octic polynomial function
  
--Position Function: 
This function will find the position vectors of all three measurements
  
  --Gibbs Functions:
This function will find the velocity that corresponds to one of the position vectors.

-Arrival_Orbit class:
Takes the input of the orbit type and orbit number from the catalog and finds the initial conditions and orbital period for the chosen orbit and sets these attributes to an object

-Transfer Function
This function is used to run the transfer design code and dependent functions. It starts by finding the manifolds and ends with the differential corrections process. Any converged solutions will be saved for the catalog and plotting. The transfer function is dependent on the following functions:
 
--Generating Manifold Initial Conditions Function:
This function will generate the initial conditions for the manifold. The use of threads will be partially helpful here since the function uses a For Loop where they are not       dependent on each other

--Processing Negative/Positive Manifolds:
These functions are used in conjunction with multithreading to speed up the process of generating the poincare section of each branch of the manifold (see Generating the Poincare Section Function)
 
--Generating the Poincare Section Function:
To ensure the manifold stops at the plane of the Earth, this function will use a basic differential corrections (single shooting) scheme to ensure the manifold will only be      propagated for the correct time. Similar to the previous function, the use of threads will be really helpful in speeding up the computational time since each branch of the manifold is independent of each other

 
--Manifold Hits the Moon Function:
This function will be used to ensure that no parts of the manifold hits the Moon. Any branch that hits the Moon will not be considered moving forward in the program
Earth-Centered Inertial (ECI) to Rotating Coordinate Frame Transformation Function:
It is assumed that the measurements of the spacecraft are found in the ECI coordinate system. Since the position and velocity will be used in later parts of the project (as a constraint and to calculate the final delta_v), it must be converted to the frame the dynamics are in, which is the rotating body frame
 
--Getting the Initial Guess for delta_v and TOF:
This function will utilize a brute force shotgunning method of finding a better initial guess for delta_v and TOF. Same as the manifold and Poincare section functions, the use of threads will be really helpful in speeding up the computational time since each combination of delta_v and TOF are independent of each other. They only need to save the final 	position of each to compare at the end of the shotgunning method. Whichever guess is closest to the position of the satellite will be used for differential corrections in the 	next function
 
--Differential Corrections Function (Inbedded in the transfer function): 
This function will be the culmination of all the previous functions. It will require the use of the position of the spacecraft as a constraint, the manifold as an initial 	position and velocity, and the guess of delta_v and TOF as initial guesses for the shooting scheme. For each iteration of the function, it will propagate the equations of 	motion, compute the error, and determine how much to change the initial conditions by. The function will end once the error is below some tolerance. The function will then rerun for the next branch of the manifold using the previous converged solution as the first initial guess. 

 
-Output Functions
These functions will serve as a way to both output and store the converged solution, output system information, and visualize the transfers

--Logger_V2:
The purpose of this function is to provide checkpoints to ensure the code is running smoothly, as well as important outputs such as for the manifold and the delta_v consuption for a converged solution. The logger will also save key pieces of information to an external log.

--Plotting_Periodic_Orbit:
The purpose of this function is to get the state of the periodic orbit after one orbital period that will be used for plotting

--Recreating_Transfer:
The purpose of this function is to recreate the transfer based on the initial conditions and delta_V from a converged solution

--Final_Outputs:
The purpose of this function is to output the catalog data to both the terminal and txt file at every converged solution. It will also output the metadata to the catalog after the first converged solutions

--Plotting_Earth_Moon_Rotating:
The purpose of this function is to transfer the Earth and Moon from the ECI to the rotating frame and plot them

--Orbit_Plot:
The purpose of this function is to plot the all the x,y, and z vectors that represent the transfers and the arrival orbit

