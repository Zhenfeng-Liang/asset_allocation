function [der] = diffDer(x, i)
% Input, x, asset price, ith direction
% Output: diff derivative matrix with respect to x_i
 
  % mean reverting model
  F = length(x);
  p = 2;               % This needs to be fixed, hard-coded right now.
  der = zeros(F, p);
  
  % lognormal model
  %lnVols = [0.33; 0.26];
  %F = length(lnVols);
  %vecD = zeros(F, 1);
  %vecD(i) = lnVols(i);
  %der = diag(vecD);
  
end

