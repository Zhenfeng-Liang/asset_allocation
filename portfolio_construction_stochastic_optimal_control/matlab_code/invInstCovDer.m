function [derx] = invInstCovDer(x, i)
% Input: x, asset price vector, with respect to ith asset, correlation matrix between underlying brownian motions
% Output: derivative of the inverse covariance matrix with respect to x_i.  i.e. (partial C_inv) / (partial x_i) 

  iCovx = invInstCov(x);
  derx = -1 * iCovx * instCovDer(x, i) * iCovx;
  
end
