classdef simprops < handle
    properties
        Ts, Tsim, Tvec, road, Vx_des, inits, inputs, ld, controlchoice, N
    end
    
    methods
        function init_conds = states_init(data1, data2,conds) 
            init_conds.Px0= conds(1); init_conds.Py0= conds(2);
            init_conds.Vx0= data1.Vx_des; init_conds.Vy0=conds(5);
            init_conds.psi0= conds(3); init_conds.psidot0= conds(6);
            init_conds.phi0= 0; init_conds.phidot0= 0;
            init_conds.theta0= 0; init_conds.thetadot0= 0;
            init_conds.omega0= data1.Vx_des/data2.Rw;
            init_conds.alpha0= 0;
        end
        
        function inputs= DeltaSineFunc(data, amplitude,frequency)
            t= (0:data.Ts:((data.Tsim-4)-data.Ts))';
            inputs = amplitude*sin(2*pi*frequency*t);
            inputs = [zeros(single(2/data.Ts),1);inputs;zeros(single(2/data.Ts),1)];
        end  
        
        function inputs= DeltaStepFunc(data,Val,StepTime)
            delta = frest.createStep('StepTime',StepTime, ...
                'StepSize',Val,'FinalTime',data.Tsim-data.Ts,'Ts',data.Ts);
            inputs=delta.Data(:,1)';
        end
        
        function tire= pacjekainit(data1)
            
                 if data1.road==1
                    % Ashpalt Front
                    tire.mu_x_f= 1.20 ;
                    tire.mu_y_f= 0.935;
                    tire.Bx_f= 11.7;
                    tire.By_f= 8.86;
                    tire.Cx_f= 1.69;
                    tire.Cy_f= 1.19;
                    tire.Ex_f= 0.377;
                    tire.Ey_f= -1.21;

                    % Ashpalt rear
                    tire.mu_x_r= 1.20;
                    tire.mu_y_r= 0.961;
                    tire.Bx_r= 11.1; 
                    tire.By_r= 9.30; 
                    tire.Cx_r= 1.69;
                    tire.Cy_r = 1.19;
                    tire.Ex_r= 0.362; 
                    tire.Ey_r= -1.11;
                
                 elseif data1.road==2
                    % Snow front
                    tire.mu_x_f= 0.407;
                    tire.mu_y_f= 0.383;
                    tire.Bx_f = 10.2; 
                    tire.By_f= 19.1;
                    tire.Cx_f= 1.96; 
                    tire.Cy_f= 0.550; 
                    tire.Ex_f= 0.651; 
                    tire.Ey_f= -2.1;

                    % Snow rear
                    tire.mu_x_r = 0.409;
                    tire.mu_y_r= 0.394;
                    tire.Bx_r= 9.71;
                    tire.By_r= 20.0;
                    tire.Cx_r= 1.96;
                    tire.Cy_r= 0.550; 
                    tire.Ex_r = 0.624;
                    tire.Ey_r= -1.93; 
                 
                 elseif data1.road==3
                    % Wet ashphalt front
                    tire.mu_x_f= 1.06;
                    tire.mu_y_f= 0.885;
                    tire.Bx_f = 12;
                    tire.By_f= 10.7;
                    tire.Cx_f= 1.80;
                    tire.Cy_f= 1.07;
                    tire.Ex_f= 0.313;
                    tire.Ey_f= -2.14;
                    
                    % Wet ashphalt rear
                    tire.mu_x_r = 1.07;
                    tire.mu_y_r= 0.911;
                    tire.Bx_r= 11.5;
                    tire.By_r= 11.3;
                    tire.Cx_r= 1.80;
                    tire.Cy_r= 1.07;
                    tire.Ex_r = 0.300;
                    tire.Ey_r= -1.97;
                 
                 else
                    % Smooth ice front
                    tire.mu_x_f= 0.172;
                    tire.mu_y_f= 0.162;
                    tire.Bx_f = 31.1;
                    tire.By_f= 28.4;
                    tire.Cx_f= 1.77;
                    tire.Cy_f= 1.48;
                    tire.Ex_f= 0.710;
                    tire.Ey_f= -1.18;
                    
                    % Smooth ice rear
                    tire.mu_x_r = 0.173;
                    tire.mu_y_r= 0.167;
                    tire.Bx_r= 29.5;
                    tire.By_r= 30;
                    tire.Cx_r= 1.77;
                    tire.Cy_r= 1.48;
                    tire.Ex_r = 0.681;
                    tire.Ey_r= -1.08;
                 end
        end
    end
end

