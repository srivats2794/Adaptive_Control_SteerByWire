classdef vehicle < handle
    properties
        g=9.81,m,Ixx,Iyy,Izz,l_f,l_r,l,w,Rw,Re,Iw,sig,h,K_theta, D_theta
        K_phif, K_phir, K_phi, D_phif,D_phir, D_phi, C_alpha_f, C_alpha_r
        epsilon, scale, linmodchoice, params
        bm %(linear bicycle model)
        tm %(pacjeka tire model)
        sr %(Steering ratio)
    end
    
    methods
        function linsys = LinearStateSpace(veh, data)
          if veh.linmodchoice==1
            A22= (-(2*(veh.C_alpha_f+veh.C_alpha_r)))/(veh.m*data.Vx_des);
            A24= -data.Vx_des-((2*((veh.C_alpha_f*veh.l_f)-(veh.C_alpha_r*veh.l_r)))/(veh.m*data.Vx_des));
            A42= (-(2*((veh.C_alpha_f*veh.l_f)-(veh.C_alpha_r*veh.l_r))))/(veh.Izz*data.Vx_des);
            A44= (-(2*((veh.l_f*veh.l_f*veh.C_alpha_f)+(veh.l_r*veh.l_r*veh.C_alpha_r))))/(veh.Izz*data.Vx_des);

            linsys.Ac= [0 1 0 0;
                        0 A22 0 A24;
                        0 0 0 1;
                        0 A42 0 A44];

            B2= (2*veh.C_alpha_f)/veh.m;
            B4= (2*veh.l_f*veh.C_alpha_f)/veh.Izz;

            linsys.Bc=[0;B2;0;B4];

            linsys.Cc= eye(4);
            linsys.Dc= zeros(4,1);

            sysc= ss(linsys.Ac,linsys.Bc,linsys.Cc,linsys.Dc);
            linsys.lin_veh= c2d(sysc,data.Ts);
            [linsys.Ad,linsys.Bd,linsys.Cd,linsys.Dd]=ssdata(linsys.lin_veh);

            linsys.om = sqrt(A22*A44-A24*A42); %natural frequency of bicycle model
            linsys.zeta= -(A22+A44)/(2*linsys.om); %damping ratio of bicycle model 
            linsys.omd = linsys.om*sqrt(1-linsys.zeta^2); %damped natural frequency of bicycle model
            linsys.Td = 2*pi/linsys.omd; %damped period of bicycle model
          
          elseif veh.linmodchoice==2
            A22= (-(2*(veh.C_alpha_f+veh.C_alpha_r)))/(veh.m*data.Vx_des);
            A23= (2*(veh.C_alpha_f+veh.C_alpha_r))/(veh.m);
            A24= (2*((-veh.C_alpha_f*veh.l_f)+(veh.C_alpha_r*veh.l_r)))/(veh.m*data.Vx_des);
            A42= (-2*((veh.C_alpha_f*veh.l_f)-(veh.C_alpha_r*veh.l_r)))/(veh.Izz*data.Vx_des);
            A43= (2*((veh.C_alpha_f*veh.l_f)-(veh.C_alpha_r*veh.l_r)))/veh.Izz;
            A44= (-2*((veh.l_f*veh.l_f*veh.C_alpha_f)+(veh.l_r*veh.l_r*veh.C_alpha_r)))/(veh.Izz*data.Vx_des);
            
            B2=  (2*veh.C_alpha_f)/veh.m;
            B4= (2*veh.l_f*veh.C_alpha_f)/veh.Izz;
            
            K21= ((2*((-veh.C_alpha_f*veh.l_f)+(veh.C_alpha_r*veh.l_r)))/(veh.m*data.Vx_des))-data.Vx_des; 
            K22= veh.g;
            K41= (-2*((veh.l_f*veh.l_f*veh.C_alpha_f)+(veh.l_r*veh.l_r*veh.C_alpha_r)))/(veh.Izz*data.Vx_des);
            
            linsys.Ac= [0 1 0 0;
                 0 A22 A23 A24;
                 0 0 0 1;
                 0 A42 A43 A44;];
            linsys.Bc= [0;B2;0;B4];
            linsys.Kc= [0 0;
                K21 K22;
                0 0;
                K41 0];
            linsys.Cc= eye(4);
            linsys.Dc= zeros(4,3);
            sysc= ss(linsys.Ac,[linsys.Bc linsys.Kc],linsys.Cc,linsys.Dc);
            linsys.lin_veh= c2d(sysc,data.Ts);
            [linsys.Ad,Bveh,linsys.Cd,linsys.Dd]=ssdata(linsys.lin_veh);
            linsys.Bd= Bveh(:,1);
            linsys.Kd_psidotdes= Bveh(:,2);
            linsys.Kd_bank= Bveh(:,3);
          
          elseif veh.linmodchoice==3
            A11= -((veh.C_alpha_f+veh.C_alpha_r)/(data.Vx_des*veh.m));
            A12= ((veh.C_alpha_r*veh.l_r-veh.C_alpha_f*veh.l_f)/(veh.m*data.Vx_des*data.Vx_des))-1;
            A21= ((-veh.C_alpha_f*veh.l_f)+(veh.C_alpha_r*veh.l_r))/(veh.Izz);
            A22= -((veh.C_alpha_f*(veh.l_f^2)+veh.C_alpha_r*(veh.l_r^2))/(veh.Izz*data.Vx_des));
            
            B1=  (veh.C_alpha_f)/(veh.m*data.Vx_des);
            B2=  (veh.l_f*veh.C_alpha_f)/veh.Izz; 
     
            linsys.Ac= [A11 A12;
                        A21 A22];
            linsys.Bc= [B1;B2];
            linsys.Cc= eye(2);
            linsys.Dc=[0;0];
            sysc= ss(linsys.Ac,linsys.Bc,linsys.Cc,linsys.Dc);
            linsys.lin_veh= c2d(sysc,data.Ts);
            [linsys.Ad,linsys.Bd,linsys.Cd,linsys.Dd]=ssdata(linsys.lin_veh); 
            
            linsys.om = sqrt(A11*A22-A12*A21); %natural frequency of bicycle model
            linsys.zeta= -(A11+A22)/(2*linsys.om); %damping ratio of bicycle model 
            linsys.omd = linsys.om*sqrt(1-linsys.zeta^2); %damped natural frequency of bicycle model
            linsys.Td = 2*pi/linsys.omd; %damped period of bicycle model
          end
        end        
    end
end
