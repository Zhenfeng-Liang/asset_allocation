function [b] = diffV(x)
% Given x, asset price, output asset dynamics diff
% Right now, I am hard coding dimension of dZ, 3, which is p dimensions in the paper
% Return: diff matrix, whose dimension is n*p
% Note: to get a dirty implementation, we let p=n at this point

  # mean reverting model
  nVols = [0.1, 0.16];        # This was hardcoded vol
  b = diag(nVols);             # b is the symbol in the paper
  
  # lognormal model
  #lnVols = [0.33; 0.26];
  #b = diag(lnVols .* x);    
  
end
