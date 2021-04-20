# Adaptive_Control_SteerByWire
MRAC (Model Reference Adaptive Control)+ Input shaping to safely control Steer-By-Wire vehicles.

Add 6DoF_plant_functions, classes and init_files folders to your Matlab Path

Run SbWAdaptiveControl file for simulation.

The timestep, road condition and operating/desired velocity can all be manipulated using the init_files/sim_params.m .

The classes are vehicle,simprops and inputshaper

vehicle class carries vehicle params and also contains functions that can initialize a desired linear model. Model 1 (using linmodchoice variable): y ydot psi psidot statespace. Model 2: ey eydot epsi epsidot statespace also known as error dynamics state space. Model 3: psidot and beta statespace also known as sideslip model.

The 6DoF_plant_functions folder contains all the files required to simulate the high DoF plant that the controller is tested against. It includes 6DoF chassis and pacjeka wheel-tire model. You can ask for next state, velocity states/states_dot and Forces as output.

Please go through the preprint.pdf file for theory. This is a confidential file and hence is not allowed to be distributed. 
