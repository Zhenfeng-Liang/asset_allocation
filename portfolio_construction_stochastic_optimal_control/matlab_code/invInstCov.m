function [invCovx] = invInstCov(x, corrMatr)
% Input: x, asset price vector
% Output: inverse matrix of the instantaneous covariance matrix

    invCovx = inv(instCov(x, corrMatr));    
end
