classdef InputShaper < handle
    properties
        SWbuffer, ismode, shaped_inputs, input_raw_history
        input_raw_vec, sys
    end
    
    methods
        function [in,input_history] = shaper(isobj, simobj,i)
            switch isobj.ismode
                case 1
                    K= exp((-isobj.sys.zeta*pi)/(sqrt(1-(isobj.sys.zeta^2))));
                    input_history = [isobj.input_raw_history, isobj.input_raw_vec(i)];
                    %first time delay @ t=0;
                    t2 = isobj.sys.Td/2; %second time delay
                    A1= 1/(1+K);
                    A2= K/(1+K);
                    input_shaped = A1*isobj.input_raw_vec(i)+A2*input_history(end-round(t2/simobj.Ts));
                    in = [isobj.shaped_inputs, input_shaped];
                    
                case 2
                    K= exp((-isobj.sys.zeta*pi)/(sqrt(1-(isobj.sys.zeta^2))));
                    input_history = [isobj.input_raw_history, isobj.input_raw_vec(i)];
                    %first time delay @ t=0;
                    t2 = isobj.sys.Td/2; %second time delay
                    t3 = isobj.sys.Td;
                    A1 = 1/(1+(2*K)+(K^2)); %amplification factor associated with first delay
                    A2 = (2*K)/(1+(2*K)+(K^2)); %amplification factor associated with second delay
                    A3= (K^2)/(1+(2*K)+(K^2));
                    input_shaped = A1*isobj.input_raw_vec(i)+A2*input_history(end-round(t2/simobj.Ts))+A3*input_history(end-round(t3/simobj.Ts));
                    in = [isobj.shaped_inputs, input_shaped];
  
                case 3
                    K= exp((-isobj.sys.zeta*pi)/(sqrt(1-(isobj.sys.zeta^2))));
                    input_history = [isobj.input_raw_history, isobj.input_raw_vec(i)];
                    %first time delay @ t=0;
                    t2 = isobj.sys.Td/2; %second time delay
                    t3 = isobj.sys.Td;
                    t4 = 1.5*isobj.sys.Td;
                    A1 = 1/(1+(3*K)+(3*(K^2))+(K^3)); %amplification factor associated with first delay
                    A2 = (3*K)/(1+(3*K)+(3*(K^2))+(K^3)); %amplification factor associated with second delay
                    A3= (3*(K^2))/(1+(3*K)+(3*(K^2))+(K^3));
                    A4= (K^3)/(1+(3*K)+(3*(K^2))+(K^3));
                    input_shaped = A1*isobj.input_raw_vec(i)+A2*input_history(end-round(t2/simobj.Ts))+A3*input_history(end-round(t3/simobj.Ts))...
                        +A4*input_history(end-round(t4/simobj.Ts));
                    in = [isobj.shaped_inputs, input_shaped];
            end
        end
    end
end