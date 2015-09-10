function [b] = diffV(x)
% Given x, asset price, output asset dynamics diff
% Right now, I am hard coding dimension of dZ, 3, which is p dimensions in the paper
% Return: diff matrix

  # mean reverting model
  nVols = [13.2, 23.0, 21.1];  # This was hardcoded vol
  b = diag(nVols);             # b is the symbol in the paper
  
  # lognormal model
  #lnVols <- c(0.33, 0.26, 0.19, 0.22, 0.31);
  #d <- diag(lnVols * x);    
  
end
