vehicle_params;
sim_params;

veh.tm = pacjekainit(data); %(pacjeka tire model)

veh.params(23)=veh.tm.mu_x_f;
veh.params(24)=veh.tm.mu_y_f;
veh.params(25)=veh.tm.Bx_f;
veh.params(26)=veh.tm.By_f;
veh.params(27)=veh.tm.Cx_f;
veh.params(28)=veh.tm.Cy_f;
veh.params(29)=veh.tm.Ex_f;
veh.params(30)=veh.tm.Ey_f;

% Ashpalt rear
veh.params(31)=veh.tm.mu_x_r;
veh.params(32)=veh.tm.mu_y_r;
veh.params(33)=veh.tm.Bx_r;
veh.params(34)=veh.tm.By_r;
veh.params(35)=veh.tm.Cx_r;
veh.params(36)=veh.tm.Cy_r;
veh.params(37)=veh.tm.Ex_r;
veh.params(38)=veh.tm.Ey_r;

veh.params(39)= veh.epsilon;
veh.params(40)= veh.scale;

Fzf= veh.m*veh.g*(veh.l_r/veh.l);
Fzr= veh.m*veh.g*(veh.l_f/veh.l);

Df= veh.tm.mu_y_f*Fzf;
Dr= veh.tm.mu_y_r*Fzr;

veh.C_alpha_f= veh.tm.By_f*veh.tm.Cy_f*Df/veh.g;
veh.C_alpha_r= veh.tm.By_r*veh.tm.Cy_r*Dr/veh.g;

clearvars Fzf Fzr Df Dr