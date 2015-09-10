function [der] = diffDer(x, i)
% Input, x, asset price, ith direction
% Output: diff derivative matrix with respect to x_i
 
  # mean reverting model
  F = 2;
  der = zeros(F, F);
  
# lognormal model
#lnVols <- c(0.33, 0.26, 0.19, 0.22, 0.31);
#F <- length(lnVols);
#vecD <- vector(mode = "double", F);
#vecD[i] <- lnVols[i];
#der <- diag(vecD);
  
end

