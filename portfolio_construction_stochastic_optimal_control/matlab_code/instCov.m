function [covx] = instCov(x, corrMatr)
% Input: x, asset price vector, and corrMatr, correlation matrix of the underlying brownian motion         
% Output: instantaneous covariance matrix of x

  dofx = diffV(x);
  covx = dofx * corrMatr * dofx';
end     

