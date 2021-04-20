function [ out ] = sign_sigmoid( in, scale )

out = -1+2/(1+exp(-scale*in));

end

