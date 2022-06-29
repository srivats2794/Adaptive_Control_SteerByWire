init.m;

 %% Uncomment for Double Lane Change Input
load('steering_ref_dlc.mat');
delta_ref= deg2rad(veh.sr*(ref_steering_angle.Data));
delta_ref= resample(delta_ref,1,(data.Ts/0.0005));
delta_ref= [delta_ref;-delta_ref;delta_ref;-delta_ref;]; % UNCOMMENT FOR REPEATED DLC
% plot(delta_ref);
clearvars ref_steering_angle
data.inputs.delta_raw= delta_ref;
[rows,~]=size(data.inputs.delta_raw);
data.N=rows;
data.Tsim= floor(data.N*data.Ts);
data.Tvec= 0:data.Ts:(data.Tsim); % Time vector

%% Setting initial conditions. 

% Passing in order Px0 Py0 Psi0 0 Vy0 Psidot0. 
% The fourth element is always 0, cuz initial Vx0 is always taken from data.Vx_des.
% Seem init_files/sim_params for context on data.Vx_des  
data.inits= states_init(data,veh,[0;0;0;0;0;0]);

%% Input Shaper Setup

veh.linmodchoice=1;
veh.bm= LinearStateSpace(veh,data);
is= InputShaper; % Declaring InputShaper class object
is.sys= LinearStateSpace(veh,data); %(bicycle model)
is.ismode= 2; % 1 for ZV, 2 for ZVD, 3 for ZVDD
is.SWbuffer = (1/data.Ts)*100;  % Steering history length
is.input_raw_history = zeros(1,is.SWbuffer); % Steering history for Time-Delay Filter

clearvars pos delta_ref rows

%% Linear reference model init
veh.linmodchoice=1;
veh.bm= LinearStateSpace(veh,data); %(bicycle model)

%% Plotter choice : 0 or 1
plot_states=1; % Do you want to plot vehicle responses?
plot_mrac_params=1; % Do you want to plot MRAC control gains?
plot_errors = 1; % Do you want to plot errors?
