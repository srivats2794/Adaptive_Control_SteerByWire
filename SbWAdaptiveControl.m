clc; clear;

IScontrol_init; % Initialization code saved in separate file. 

%% MRAC control params setup
Kdel(1)=0; Ky(:,1)=zeros(2,1);Ke(:,1)= zeros(2,1);

gamma_y=diag([0.001;0.1]); gamma_del=10;  gamma_e= 4*eye(2);%1*ones(2,1);
Q_mrac= eye(2);
A= [veh.bm.Ac(2,2) veh.bm.Ac(2,4);
     veh.bm.Ac(4,2) veh.bm.Ac(4,4)];
B= [veh.bm.Bc(2);veh.bm.Bc(4)];
C=eye(2);
D= zeros(2,1);
sys_mrac= ss(A,B,C,D);
sysd_mrac= c2d(sys_mrac,data.Ts);
[Am,Bm,~,~]=ssdata(sysd_mrac);
P_mrac= dlyap(Am,Q_mrac);
B_mrac= 4e-3*[1;1];

%% Vehicle response recorder 
% Pre-Allocation to speed up code. 18 states of plant. is subscript is for inputshaper only
% mrac subscript is for mrac control. ref subscript is for reference model.
states_history_is= zeros(18,data.N+1); 
states_history_mrac= zeros(18,data.N+1);
states_history_ref= zeros(2,data.N+1); % ydot psidot of reference model
pos_history_ref= zeros(3,data.N+1); % x y and psi of reference model

%% Initial state vectors
states_plant_is= [data.inits.Px0;data.inits.Py0; data.inits.psi0; data.inits.Vx0;data.inits.Vy0; ...
         data.inits.psidot0; data.inits.theta0; data.inits.thetadot0; data.inits.phi0; ...
         data.inits.phidot0; data.inits.omega0; data.inits.omega0; data.inits.omega0; ...
         data.inits.omega0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0];
vxy_plant_is=[data.inits.Vx0;data.inits.Vy0]; 
states_history_is(:,1)=states_plant_is;

states_plant_mrac= states_plant_is;
vxy_plant_mrac=vxy_plant_is; 
states_history_mrac(:,1)=states_plant_mrac;

states_ref = [data.inits.Vy0;data.inits.psidot0];
pos_history_ref(:,1)= [data.inits.Px0;data.inits.Py0;data.inits.psi0];
states_history_ref(:,1)= [data.inits.Vy0;data.inits.psidot0];

meas=[data.inits.Vy0;data.inits.psidot0];
vxy_ref=[data.inits.Vx0;data.inits.Vy0];
err=[0;0];

for i= 1:data.N
    is.input_raw_vec(i)= data.inputs.delta_raw(i);
    [is.shaped_inputs,is.input_raw_history]= shaper(is,data,i); % Function to find shaped inputs
      
    %% Final inputs recorded ready to be used with dynamics propagation
    data.inputs.delta_ref(i) = is.shaped_inputs(i);
    data.inputs.delta_is(i) = is.shaped_inputs(i); % data.inputs.theta_sw(i);
    data.inputs.delta_mrac(i)= is.shaped_inputs(i)*Kdel(i)+ ...
        Ky(:,i)'*meas+Ke(:,i)'*err;
 
    %% Reference Model (Linear plant) states propagation
    pos_history_ref(:,i+1)= pos_history_ref(:,i)+data.Ts*[vxy_ref;states_ref(2)]; 
    states_ref = Am*states_ref+Bm*data.inputs.delta_ref(i);
    %yddot_lin= states_ref_dot(2)+data.Vx_des*states_ref(4);
    rot=[cos(pos_history_ref(3,i+1)) -sin(pos_history_ref(3,i+1));
         sin(pos_history_ref(3,i+1))  cos(pos_history_ref(3,i+1))];
    vxy_ref= rot*[data.Vx_des; states_ref(1)];
    states_history_ref(2,i+1)= states_ref(2);
    states_history_ref(1,i+1)= vxy_ref(2);
        
    %% Nonlinear plant forward propagation   
    
    ins_is= [0;0;data.inputs.delta_is(i)];
    ins_mrac= [0;0;data.inputs.delta_mrac(i)];
        
    [states_plant_is,vel_vec_is] = plant_euler_forward(data.Ts,states_plant_is,ins_is,veh.params);
    states_plant_mrac = plant_euler_forward(data.Ts,states_plant_mrac,ins_mrac,veh.params);
    
    % Rotation matrices to rotate velocities from dynamics propagation from
    % body to world
    rotP1= [cos(states_plant_is(3)) -sin(states_plant_is(3));
           sin(states_plant_is(3))  cos(states_plant_is(3))];
    
    rotP2= [cos(states_plant_mrac(3)) -sin(states_plant_mrac(3));
            sin(states_plant_mrac(3))  cos(states_plant_mrac(3))];
    
    vxy_plant_is=rotP1*[states_plant_is(4);states_plant_is(5)];
    vxy_plant_mrac=rotP2*[states_plant_mrac(4);states_plant_mrac(5)];
    
    %% Recording states in the history vector for plotting purposes
    states_history_is(:,i+1)= states_plant_is;
    states_history_mrac(:,i+1)= states_plant_mrac;
        
    states_history_is(4:5,i+1)= vxy_plant_is;
    states_history_mrac(4:5,i+1)= vxy_plant_mrac;
    
    %% MRAC params
    ydot_ref= vxy_ref(2);
    ydot_plant= states_history_mrac(5,i+1);
        
    % MRAC Errors
    e_ydot= -(ydot_ref-ydot_plant);
    e_psidot= -(states_ref(2)-states_plant_mrac(6));
    err=[e_ydot;e_psidot];
    mrac_err(:,i)=err;
    
    meas=[ydot_ref;states_ref(2)];
    
    % Gains
    Ky_dot= -gamma_y*meas*err'*P_mrac*B_mrac;
    Kdel_dot= -gamma_del*data.inputs.delta_ref(i)*err'*P_mrac*B_mrac;
    Ke_dot= -gamma_e*(err*err')*P_mrac*B_mrac; %-(err*err')*gamma_e;
    
    Kdel(i+1)= Kdel_dot*data.Ts+Kdel(i);
    Ky(:,i+1)= Ky_dot*data.Ts+Ky(:,i);
    Ke(:,i+1)=Ke_dot*data.Ts+Ke(:,i);
end

beta_lqr= atan2(states_history_is(5,:),states_history_is(4,:));

states_history_is(19,:)= states_history_is(3,:)-beta_lqr;

beta_mrac= atan2(states_history_mrac(5,:),states_history_mrac(4,:));

states_history_mrac(19,:)= states_history_mrac(3,:)-beta_mrac;

clearvars -except states_history_is states_history_mrac ref_states data veh ...
                  pos_history_ref states_history_ref plot_states ...
                  plot_mrac_params K_psi Kdel Ky mrac_err Ke plot_errors

[~,columns]=size(states_history_is);
Tvec_sim= 0:data.Ts:((columns-1)*data.Ts); % Time vector

if plot_states==1
    figure(11)
    nexttile
    plot(states_history_is(1,:),states_history_is(2,:),'-k', ...
        pos_history_ref(1,:),pos_history_ref(2,:),'--r', ...
        states_history_mrac(1,:),states_history_mrac(2,:),':b','LineWidth',4);
    set(gca,'fontsize',20);
    title ('Trajectory','FontSize',25);
    legend('Input Shaper','Reference','MRAC','FontSize',20);
    ylabel('Lateral Displacement(m)','FontSize',25); xlabel('Longitudinal Displacement(m)','FontSize',25);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_is(3,:)),'-k', ...
         Tvec_sim,rad2deg(pos_history_ref(3,:)),'--r', ...
         Tvec_sim,rad2deg(states_history_mrac(3,:)),':c','LineWidth',2);
    set(gca,'fontsize',12);
    title ('Yaw','FontSize',10);
    legend('IS','Linear Plant','MRAC','FontSize',7,'Location','southwest');
    ylabel('$\psi$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_is(6,:)),'-k', ...
         Tvec_sim,rad2deg(states_history_ref(2,:)),'--r', ...
         Tvec_sim,rad2deg(states_history_mrac(6,:)),':b','LineWidth',4);
    set(gca,'fontsize',20);
    title ('Yaw Rate','FontSize',25);
    legend('IS','Linear Plant','MRAC','FontSize',7,'Location','southwest');
    ylabel('$\dot{\psi}$ (deg)','FontSize',30,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',25);
    
    nexttile
    plot(Tvec_sim,states_history_is(5,:),'-k', ...
        Tvec_sim,states_history_ref(1,:),'--r', ...
        Tvec_sim,states_history_mrac(5,:),':c','LineWidth',2);
    set(gca,'fontsize',12);
    title ('Lateral Velocity','FontSize',10);
    legend('IS','Linear Plant','MRAC','FontSize',7,'Location','southwest');
    ylabel('$\dot{y}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',25);
    
    nexttile
    plot(Tvec_sim,states_history_is(4,:),'-k', ...
        Tvec_sim,data.Vx_des*ones(1,columns),'--r', ...
        Tvec_sim,states_history_mrac(4,:),':c','LineWidth',2);
    set(gca,'fontsize',12);
    title ('Longitudinal Velocity','FontSize',10);
    legend('IS','Linear Plant','MRAC','FontSize',7,'Location','southwest');
    ylabel('$\dot{x}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_is(19,:)),'-k', ...
         Tvec_sim,rad2deg(states_history_mrac(19,:)),':c','LineWidth',2);
    set(gca,'fontsize',12);
    title ('SideSlip','FontSize',10);
    legend('IS','MRAC','FontSize',7,'Location','southwest');
    ylabel('$\beta$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
        
    sgtitle('VEHICLE RESPONSE PLOTS','FontSize',20)

    figure(32)
    plot(Tvec_sim(1:end-1),rad2deg(data.inputs.delta_is),'-k', ...
        Tvec_sim(1:end-1),rad2deg(data.inputs.delta_raw),'--r', ...
        Tvec_sim(1:end-1),rad2deg(data.inputs.delta_mrac),':b', 'LineWidth',4);
    set(gca,'fontsize',20);
    %xlim([0 101])
    title ('Control inputs','FontSize',25);
    legend('Shaped Input','Raw Input','MRAC','FontSize',20,'Location','southwest');
    ylabel('$\delta$ (deg)','FontSize',30,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',25);   
end

if plot_mrac_params==1
    figure(13)
    nexttile
    plot(Tvec_sim(1:end-1),mrac_err(1,:),'-k', 'LineWidth',2);
    title ('Lateral Velocity error','FontSize',15);
    ylabel('$\varepsilon_{\dot{y}}$ (m/s)','FontSize',20,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    xlim([0 101])
    ylim([-1 1])
    nexttile
    plot(Tvec_sim(1:end-1),mrac_err(2,:),'-k', 'LineWidth',2);
    title ('Yaw Rate error','FontSize',15);
    ylabel('$\varepsilon_{\dot{\psi}}$ (rad/s)','FontSize',20,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    xlim([0 101])
    ylim([-1 1])
end

if plot_errors==1
    figure(14)
    nexttile
    plot(Tvec_sim,Ky(1,:),'-k', 'LineWidth',2);
    title('','FontSize',10);
    nexttile
    plot(Tvec_sim,Ky(2,:),'-k', 'LineWidth',2);
    title('Psidot state gain','FontSize',10);
    nexttile
    plot(Tvec_sim,Ke(1,:),'-k', 'LineWidth',2);
    title('Ydot error gain','FontSize',10);
    nexttile
    plot(Tvec_sim,Ke(2,:),'-k', 'LineWidth',2);
    title('Psidot error gain','FontSize',10);
    nexttile
    plot(Tvec_sim,Kdel,'-k', 'LineWidth',2);
    title('Input feedforward gain','FontSize',10);
end