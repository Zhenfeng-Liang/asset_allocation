function [derx] = instCovDer2(x, i, j, corrMatr)
% Input: x, asset price vector, with respect to ith and jth asset
% Output: second derivative of covariance matrix with respect to ith and jth asset

  dofx = diffV(x);
  ddofxi = diffDer(x, i);
  ddofxj = diffDer(x, j);
  dd2ofxij = diffDer2(x, i, j);
 
  derx = dd2ofxij * corrMatr * dofx' + ddofxi * corrMatr * ddofxj' + ddofxj * corrMatr * ddofxi' + dofx * corrMatr * diffDer2(x, i, j)';
end
