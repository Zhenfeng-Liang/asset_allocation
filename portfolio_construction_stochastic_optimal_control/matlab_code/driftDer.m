function [der] = driftDer(x, i)
% Input: x, asset price vector, with respect to ith asset price
% Output: drift vector derivative with respect to x_i
 
  # mean reverting model
  lambda = [21.0, 13.2];
  F = length(lambda);
  der = zeros(1, F);
  der(i) = -lambda(i);
  
  # lognormal model
  #mu <- c(0.4, 0.12, 0.05, 0.16, 0.02);
  #F <- length(mu);
  #der <- vector(mode = "double", F);
  #der[i] <- mu[i];
  
end

