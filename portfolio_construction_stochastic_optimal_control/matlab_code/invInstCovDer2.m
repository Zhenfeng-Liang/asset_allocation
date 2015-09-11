function [derx] = invInstCovDer2(x, i, j, corrMatr)
% Input: x, asset price vector, with respect to ith and jth assets
% Output: second derivative of inverse covariance matrix with respect to ith and jth asset

  iCovx = invInstCov(x, corrMatr);
  covxdi = instCovDer(x, i, corrMatr);
  covxdj = instCovDer(x, j, corrMatr);
  derx = iCovx * (covxdi * iCovx * covxdj + covxdj * iCovx * covxdi - instCovDer2(x, i, j, corrMatr)) * iCovx;
  
end
