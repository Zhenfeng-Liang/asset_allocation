function [covx] = instCov(x)
% Input: x, asset price vector, and corrMatr, correlation matrix of the underlying brownian motion         
% Output: instantaneous covariance matrix of x
  
  global corrMatr;
  dofx = diffV(x);
  covx = dofx * corrMatr * dofx';
end     

