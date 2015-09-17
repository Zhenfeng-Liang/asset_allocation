function [der] = driftDer(x, i)
% Input: x, asset price vector, with respect to ith asset price
% Output: drift vector derivative with respect to x_i
 
  # mean reverting model
  lambda = [21.0, 13.2];
  F = length(lambda);
  der = zeros(F, 1);
  der(i) = -lambda(i);
  
  # lognormal model
  #mu = [0.4, 0.12];
  #F = length(mu);
  #der = zeros(F, 1);
  #der(i) <- mu(i);
  
end

