function [derx] = instCovDer(x, i, corrMatr)
% Input: x, asset price vector, derivative with respect to ith asset, corrMatr, correlation matrix between underlying Brownian Motion
% Output: derivative of covariance matrix

  dofx = diffV(x);
  ddofx = diffDer(x, i);
  derx = ddofx * corrMatr * dofx' + dofx * corrMatr * ddofx';

end
