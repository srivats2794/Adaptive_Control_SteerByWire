veh=vehicle;

veh.sr= 1/20; veh.params(1)=veh.sr;% Steering ratio

veh.m = 2100;   veh.params(2)=veh.m;% mass of car

veh.Ixx = 765; veh.params(3)=veh.Ixx;

veh.Iyy = 3477; veh.params(4)=veh.Iyy;

veh.Izz = 3900;veh.params(5)=veh.Izz;   % Moment of inertia
veh.l_f = 1.3; veh.params(6)=veh.l_f;  % Front axle to center of gravity
veh.l_r = 1.5; veh.params(7)=veh.l_r; % Rear axle to center of gravity
veh.l= veh.l_f+veh.l_r; veh.params(8)=veh.l;
veh.w = 0.8;  veh.params(9)=veh.w;  % width (from wheel to center)
veh.Rw = 0.3; veh.Re=veh.Rw; veh.params(10)=veh.Rw;% Tire effective radius
veh.Iw = 4.0; veh.params(11)=veh.Iw;  % Tire spinning inertia
veh.sig = 0.3; veh.params(12)=veh.sig; % relaxation length
veh.h = 0.5;   veh.params(13)=veh.h; % distance from roll to mass center
% 
veh.K_theta = 363540; veh.params(14)=veh.K_theta;
veh.D_theta = 30960; veh.params(15)=veh.D_theta;
% 
veh.K_phif = 89000; veh.params(16)=veh.K_phif;
veh.K_phir = 89000; veh.params(17)=veh.K_phir;
veh.K_phi = veh.K_phif+veh.K_phir; veh.params(18)=veh.K_phi;
veh.D_phif = 8000; veh.params(19)=veh.D_phif;
veh.D_phir = 8000; veh.params(20)=veh.D_phir;
veh.D_phi = veh.D_phif+veh.D_phir; veh.params(21)=veh.D_phi;
veh.params(22)= veh.g;
% Cornering stiffness initial estimate given by
% 16-17% of the load on the tire. 
% 
% veh.C_alpha_f= veh.m*0.165*(veh.l_f/veh.l)*57.2958;
% veh.C_alpha_r= veh.m*0.165*(veh.l_r/veh.l)*57.2958;

veh.scale = 1e3;
veh.epsilon = 1e-6;

