function [val] = hamDxx(x, p)
% Input: x, asset price vector, p momentum vector
% Output: a F x F matrix, second derivative of hamilton with respect to x_ij

  global oneOverGamma;
  global kappa;

  a = driftV(x);  
  
  F = length(x);
  term1 = zeros(F, F);
  term2 = zeros(F, F);
  term3 = zeros(F, F);
  
  for i = 1:F      
    for j = 1:F       
      term1(i,j) = 0.5 * oneOverGamma * p' * instCovDer2(x, i, j) * p;
      term2(i,j) = oneOverGamma * p' * driftDer2(x, i, j);
      term3(i,j) = kappa * ( a' * invInstCov(x) * driftDer2(x, i, j) ...
                            + driftDer(x, j)' * invInstCov(x) * driftDer(x, i) ...
                            + a' * invInstCovDer(x, j) * driftDer(x, i) ...
                            + a' * invInstCovDer(x, i) * driftDer(x, j) ...
                            + 0.5 * a' * invInstCovDer2(x, i, j) * a);
      
      val = term1 + term2 + term3;
    end
  end
end

