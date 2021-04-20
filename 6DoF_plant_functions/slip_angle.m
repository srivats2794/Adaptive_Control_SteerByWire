function [ alpha ] = slip_angle( vl, vc, epsilon, scale )

alpha = atan(sign_sigmoid(vl,scale)*vc/sqrt(vl.^2+epsilon));
% alpha = atan(vc/vl);

end

