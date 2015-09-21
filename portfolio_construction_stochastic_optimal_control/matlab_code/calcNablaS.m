function [nablaS] = calcNablaS(x, t, T, timeStep, tol, maxIter)
% Input: x: initial asset price vector, t: starting time, T: terminal time, timeStep: time interval between points, tol: tolerance, maxIter: maximum iteration times.
% Output: gradient of S with respect to x vector
% Note: make sure (T - t) / timeStep is odd. 

  F = length(x);
  S = calcS(x, t, T, timeStep, tol, maxIter);
  nablaS = zeros(F, 1);
  eps = 1e-5;

  for i = 1:F
      xEps = x;
      xEps(i) = xEps(i) + eps;
      nablaS(i) = (calcS(xEps, t, T, timeStep, tol, maxIter) - S) / eps;
  end

end
