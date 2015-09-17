function [der2] = diffDer2(x, i, j)
% Input, x, asset price vector, with respect to ith and jth asset derivatives
% Output, second derivative of diff matrix with respect to ith and jth asset

  % mean reverting and lognormal models
  F = length(x);
  p = length(x); % This needs to be fixed. Normally, this is not the same as n. To simplify the problem, let it equal at this point.

  der2 = zeros(F, p);   
  
end
