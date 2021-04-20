function [ slip ] = slip_ratio( r, v, omega, epsilon, scale )

slip = sign_sigmoid(v,scale)*(r.*omega-v)./(sqrt(v.^2+epsilon));
% slip = (r.*omega-v)./(v);

end

