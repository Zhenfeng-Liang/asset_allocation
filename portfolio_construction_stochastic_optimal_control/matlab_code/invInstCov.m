function [invCovx] = invInstCov(x)
% Input: x, asset price vector
% Output: inverse matrix of the instantaneous covariance matrix
    
    invCovx = inv(instCov(x));    
end
