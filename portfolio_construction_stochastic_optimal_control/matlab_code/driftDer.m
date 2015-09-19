function [der] = driftDer(x, i)
% Input: x, asset price vector, with respect to ith asset price
% Output: drift vector derivative with respect to x_i
 
  global isMeanReverting;

  if isMeanReverting

    % mean reverting model
    global lambda;
    F = length(lambda);
    der = zeros(F, 1);
    der(i) = -lambda(i);

  else
    % lognormal model
    global lnMu;
    F = length(lnMu);
    der = zeros(F, 1);
    der(i) <- lnMu(i);
  end

end

