function [a] = driftV(x)
% Given x, asset price, output asset dynamics drift
% Right now, I am hard coding two dimensions drift and lambda
% Return: drift vector

  global isMeanReverting;

  if isMeanReverting
    % mean reverting model
    global MRMu;
    global lambda;
    a = lambda .*  (MRMu - x);  % a is the symbol on paper
  else
    % lognormal model
    global lnMu;
    a = lnMu .* x;
  end                     
  
end
