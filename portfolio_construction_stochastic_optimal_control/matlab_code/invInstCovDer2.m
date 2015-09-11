function [derx] = invInstCovDer2(x, i, j)
% Input: x, asset price vector, with respect to ith and jth assets
% Output: second derivative of inverse covariance matrix with respect to ith and jth asset

  iCovx = invInstCov(x);
  covxdi = instCovDer(x, i);
  covxdj = instCovDer(x, j);
  derx = iCovx * (covxdi * iCovx * covxdj + covxdj * iCovx * covxdi - instCovDer2(x, i, j)) * iCovx;
  
end
