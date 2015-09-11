function [val] = hamDx(x, p)
% Input: x, asset price vector, p, momentum vector
% Output: gradient of hamilton with respect to x

    global oneOverGamma;
    global kappa;

    a = driftV(x);  
    F = length(x);

    term1 = zeros(1, F);
    term2 = zeros(1, F);
    term2 = zeros(1, F);
    
    for i = 1:F    
      term1(i) = 0.5 * oneOverGamma * p' * instCovDer(x, i) * p;
      term2(i) = oneOverGamma * p' * driftDer(x, i);
      term3(i) = kappa * a' * (invInstCov(x) * driftDer(x, i) + 0.5 * invInstCovDer(x, i) * a);  
    end
    
    val = term1 + term2 + term3;
end
