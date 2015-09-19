function [b] = diffV(x)
% Given x, asset price, output asset dynamics diff
% Return: diff matrix, whose dimension is n*p
% Note: to get a dirty implementation, we let each asset price only
% have one Brownian Motion. So b is generated as a diagonal matrix.

  global isMeanReverting

  if isMeanReverting
    
    % mean reverting model
    global nVols;
    b = diag(nVols);             % b is the symbol in the paper
  
  else    
  
    % lognormal model
    global lnVols;
    b = diag(lnVols .* x);    
  
  end
end
