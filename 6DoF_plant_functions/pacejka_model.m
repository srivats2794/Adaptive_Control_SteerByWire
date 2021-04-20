function [Fl,Fc] = pacejka_model(alpha, s, Fz, front, epsilon, pacejka_params)  
   s0 = s;
   alpha0 = alpha;
   
    mu_x_f = pacejka_params(1);
    Bx_f = pacejka_params(3);
    Cx_f = pacejka_params(5);
    Ex_f = pacejka_params(7);
    
    mu_y_f = pacejka_params(2);
    By_f = pacejka_params(4);
    Cy_f = pacejka_params(6);
    Ey_f = pacejka_params(8);
    
    mu_x_r = pacejka_params(9);
    Bx_r = pacejka_params(11);
    Cx_r = pacejka_params(13);
    Ex_r = pacejka_params(15);
    
    mu_y_r = pacejka_params(10);
    By_r = pacejka_params(12);
    Cy_r = pacejka_params(14);
    Ey_r = pacejka_params(16);
   
   if front == 1
       mu_x = mu_x_f;
       Bx = Bx_f;
       Cx = Cx_f;
       Ex = Ex_f;
       
       mu_y = mu_y_f;
       By = By_f;
       Cy = Cy_f;
       Ey = Ey_f;
   else
       mu_x = mu_x_r;
       Bx = Bx_r;
       Cx = Cx_r;
       Ex = Ex_r;
       
       mu_y = mu_y_r;
       By = By_r;
       Cy = Cy_r;
       Ey = Ey_r;
   end

   % NOMINAL Longitudinal tire force
   Fx0 = mu_x*Fz*sin(Cx*atan(Bx*(1-Ex)*s0+Ex*atan(Bx*s0)));
   
   % NOMINAL Lateral tire force
   Fy0 = mu_y*Fz*sin(Cy*atan(By*(1-Ey)*alpha0+Ey*atan(By*alpha0)));
   
   Fl = Fx0;

   Fc = Fy0*sqrt(1-(Fx0/(mu_x*Fz))^2+epsilon);
   
end
