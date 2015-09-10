function [der] = diffDer(x, i)
% Input, x, asset price, ith direction
% Output: diff derivative matrix with respect to x_i
 
  # mean reverting model
  F = length(x);
  p = 2;               # This needs to be fixed, hard-coded right now.
  der = zeros(F, p);
  
# lognormal model
#lnVols <- c(0.33, 0.26, 0.19, 0.22, 0.31);
#F <- length(lnVols);
#vecD <- vector(mode = "double", F);
#vecD[i] <- lnVols[i];
#der <- diag(vecD);
  
end

