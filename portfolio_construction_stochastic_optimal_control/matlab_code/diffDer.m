function [der] = diffDer(x, i)
% Input, x, asset price, ith direction
% Output: diff derivative matrix with respect to x_i
 
  global isMeanReverting;
  
  if isMeanReverting
    
    % mean reverting model
    F = length(x);
    p = size(diffV(x), 2);      % p is the number of column of diff
    der = zeros(F, p);

  else
    
    % lognormal model
    global lnVols;
    F = length(lnVols);
    vecD = zeros(F, 1);
    vecD(i) = lnVols(i);
    der = diag(vecD);
  
  end
end

