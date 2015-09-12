function [val] = hamDxp(x, p)
% Input: x, asset price vector, p, momentum vector
% Output: second derivative matrix of Hamilton with respect to x and p
% Note: ith row of the result matrix is the first derivative with respect to x_i, jth column is the first derivative with respect to p_j

  global oneOverGamma;
  F = length(x);
  val = zeros(F, F);
  
  for i = 1:F  
     val(i,:) = oneOverGamma * (instCovDer(x, i) * p + driftDer(x, i));
  end
end

