function [xnext,vel_curr,forces]= plant_euler_forward(t,x,u,p)
    Neul= t/0.001;
    for i=1:Neul
    % 18 States (in order) - Px,Py,Psi,Vx,Vy,PsiDot,Theta,ThetaDot, Phi, ...
    % PhiDot, omega1,2,3,4, alpha1,2,3,4 
    % 1-LF, 2-RF, 3-LR, 4-RR
      psi= x(3);
      Vx= x(4); Vy= x(5); psidot= x(6);
      theta= x(7); thetadot= x(8);
      phi= x(9); phidot= x(10);
      omega1= x(11); omega2= x(12);
      omega3= x(13); omega4= x(14);
      alpha1= x(15); alpha2= x(16);
      alpha3= x(17); alpha4= x(18);
      Tf= u(1); Tr= u(2); delta= u(3);
      
      %% Solving Fz (4 eqns 4 unknowns)
      tau_theta= (p(14)*theta)+(p(15)*thetadot);
      tau_phi= (p(18)*phi)+(p(21)*phidot);
      tau_phif= (p(16)*phi)+(p(19)*phidot);
      tau_phir= (p(17)*phi)+(p(20)*phidot);
      
      Fzf = (p(7)*p(2)*p(22) + tau_theta)/(p(8));
      Fzr = p(2)*p(22) - Fzf;
      
      Fz1 = 1/2*(Fzf - 1/p(9)*(tau_phif));
      Fz2 = Fzf - Fz1;
      Fz3 = 1/2*(Fzr - 1/p(9)*(tau_phir));
      Fz4 = Fzr - Fz3;
      %% Vx and Vy at each wheel
      rot1= [cos(delta) -sin(delta);
             sin(delta) cos(delta)];
      Vx1_ur= Vx-(p(9)*psidot); Vx2_ur= Vx+(p(9)*psidot);
      Vx3_ur= Vx1_ur; Vx4_ur= Vx2_ur;
      Vy1_ur= Vy+(p(6)*psidot); Vy3_ur= Vy-(p(7)*psidot);
      Vy2_ur= Vy1_ur; Vy4_ur= Vy3_ur;
      Vt1 = rot1*[Vy1_ur; Vx1_ur];
      Vy1= Vt1(1); Vx1= Vt1(2);
      Vt2 = rot1*[Vy2_ur; Vx2_ur];
      Vy2= Vt2(1); Vx2= Vt2(2);
      Vy3= Vy3_ur; Vx3= Vx3_ur;
      Vy4= Vy4_ur; Vx4= Vx4_ur;
      
      %% Slip Angles
      alpha1_aux= slip_angle(Vx1,Vy1,p(39),p(40)); %(-atan(Vy1/Vx1)-alpha1)*(Vx1/p(12));
      alpha2_aux= slip_angle(Vx2,Vy2,p(39),p(40));%(-atan(Vy2/Vx2)-alpha2)*(Vx2/p(12));
      alpha3_aux= slip_angle(Vx3,Vy3,p(39),p(40));%(-atan(Vy3/Vx3)-alpha3)*(Vx3/p(12));
      alpha4_aux= slip_angle(Vx4,Vy4,p(39),p(40));%(-atan(Vy4/Vx4)-alpha4)*(Vx4/p(12));
      
      alphadot1= -Vx1/p(12)*(alpha1+alpha1_aux);
      alphadot2= -Vx2/p(12)*(alpha2+alpha2_aux);
      alphadot3= -Vx3/p(12)*(alpha3+alpha3_aux);
      alphadot4= -Vx4/p(12)*(alpha4+alpha4_aux);
      
      %% Slip Ratios
      kappa1= slip_ratio(p(10),Vx1,omega1,p(39),p(40));%((p(10)*omega1)-Vx1)/Vx1;
      kappa2= slip_ratio(p(10),Vx2,omega2,p(39),p(40));%((p(10)*omega2)-Vx2)/Vx2;
      kappa3= slip_ratio(p(10),Vx3,omega3,p(39),p(40));%((p(10)*omega3)-Vx3)/Vx3;
      kappa4= slip_ratio(p(10),Vx4,omega4,p(39),p(40));%((p(10)*omega4)-Vx4)/Vx4;
      
      %% Magic Formula
      pacejka_params= p(23:38);
      [Fx1_ur, Fy1_ur] = pacejka_model(alpha1, kappa1, Fz1, 1, p(39), pacejka_params);  % front
      [Fx2_ur, Fy2_ur] = pacejka_model(alpha2, kappa2, Fz2, 1, p(39), pacejka_params);  % front
      [Fx3, Fy3] =       pacejka_model(alpha3, kappa3, Fz3, 0, p(39), pacejka_params);  % rear
      [Fx4, Fy4] =       pacejka_model(alpha4, kappa4, Fz4, 0, p(39), pacejka_params);  % rear
      
      Fxy1 = rot1*[Fx1_ur; Fy1_ur];
      Fx1 = Fxy1(1); Fy1 = Fxy1(2);
      Fxy2 = rot1*[Fx2_ur; Fy2_ur];
      Fx2 = Fxy2(1); Fy2 = Fxy2(2);
      
      %% Omega Dots
      omegadot1= (Tf/2-(p(10)*Fx1))/p(11);
      omegadot2= (Tf/2-(p(10)*Fx2))/p(11);
      omegadot3= (Tr/2-(p(10)*Fx3))/p(11);
      omegadot4= (Tr/2-(p(10)*Fx4))/p(11);
      
      %% Net Forces and Moments
      Fx= Fx1+Fx2+Fx3+Fx4;
      Fy= Fy1+Fy2+Fy3+Fy4;
      Mz= p(6)*(Fy1+Fy2) + p(9)*(Fx2-Fx1) - p(7)*(Fy3+Fy4) - p(9)*(Fx4+Fx3);   
      
      %% Euler Angles' Angular Accelerations
      
      num1= -(tau_phi) + p(13)*(Fy*cos(phi)*cos(theta)+p(2)*p(22)*sin(phi)) + ...
          psidot*(p(4)-p(5))*(psidot*sin(phi)*cos(phi)*cos(theta)+ ...
          phidot*sin(theta)*sin(phi)*cos(phi)) + psidot*thetadot*(cos(phi)^2*p(4)+sin(phi)^2*p(5));
      
      den1= (p(3)*cos(theta)^2+p(4)*sin(theta)^2*sin(phi)^2+p(5)*sin(theta)^2*cos(phi)^2);
      phiddot= num1/den1;
      
      num2= -(tau_theta) + p(13)*(p(2)*p(22)*sin(theta)*cos(phi)-Fx*cos(theta)*cos(phi)) + ...
          psidot*(psidot*sin(theta)*cos(theta)*(p(3)-p(4)+cos(phi)^2*(p(4)-p(5)))- ...
          phidot*(cos(theta)^2*p(3)+sin(phi)^2*sin(theta)^2*p(4)+sin(theta)^2*cos(phi)^2*p(5))- ...
          thetadot*(sin(theta)*sin(phi)*cos(phi)*(p(4)-p(5))));
      den2= (p(4)*cos(phi)^2+p(5)*sin(phi)^2);
      thetaddot= num2/den2;
      
      num3= Mz - p(13)*(Fx*sin(phi)+Fy*sin(theta)*cos(phi));
      den3= (p(3)*sin(theta)^2+cos(theta)^2*(p(4)*sin(phi)^2+p(5)*cos(phi)^2));
      
      psiddot= num3/den3;
      %% Net Linear Accelerations
      Vxdot= Vy*psidot + p(13)*(sin(theta)*cos(phi)*(psidot^2+phidot^2+thetadot^2)- ...
          sin(phi)*psiddot-2*cos(phi)*phidot*psidot-cos(theta)*cos(phi)*thetaddot+ ...
          2*cos(theta)*sin(phi)*thetadot*phidot+sin(theta)*sin(phi)*phiddot)+Fx/p(2);
      
      Vydot= -Vx*psidot + p(13)*(-sin(theta)*cos(phi)*psiddot-sin(phi)*psidot^2- ...
          2*cos(theta)*cos(phi)*thetadot*psidot+ ...
          sin(theta)*sin(phi)*phidot*psidot-sin(phi)*phidot^2+cos(phi)*phiddot) + Fy/p(2);
      
      %% Forward propagation
      % 18 States (in order) - Px,Py,Psi,Vx,Vy,PsiDot,Theta,ThetaDot, Phi, ...
      % PhiDot, omega1,2,3,4, alpha1,2,3,4
      rot= [cos(psi) -sin(psi);
          sin(psi)  cos(psi)];
      Vxy= rot*[Vx;Vy];
      Vx= Vxy(1); Vy= Vxy(2);
      dx= [Vx;Vy;psidot;Vxdot;Vydot;psiddot;thetadot;thetaddot;phidot;phiddot; ...
          omegadot1;omegadot2;omegadot3;omegadot4;alphadot1;alphadot2;alphadot3;alphadot4];
      x= dx*0.001+x;
    end
    
    if nargout == 2
       vel_curr= dx;
    elseif nargout == 3
       vel_curr= dx;
       forces.Fx= [Fx1;Fx2;Fx3;Fx4];
       forces.Fy= [Fy1;Fy2;Fy3;Fy4];
       forces.Fz= [Fz1;Fz2;Fz3;Fz4];
    end
    xnext=x;
end