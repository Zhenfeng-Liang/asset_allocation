function [val] = hamDp(x, p)
% Input: x, asset price vector, p, momentum vector
% Output: derivative of Hamilton with respect to momentum

    global oneOverGamma;
    val = oneOverGamma * (instCov(x) * p + driftV(x));

end 
