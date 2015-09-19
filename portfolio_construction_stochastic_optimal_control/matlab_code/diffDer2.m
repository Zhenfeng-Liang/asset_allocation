function [der2] = diffDer2(x, i, j)
% Input, x, asset price vector, with respect to ith and jth asset derivatives
% Output, second derivative of diff matrix with respect to ith and jth asset

  % mean reverting and lognormal models
  F = length(x);
  p = size(diffV(x), 2);
  der2 = zeros(F, p);   
  
end
