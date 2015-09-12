function [val] = hamDpp(x, p)
% Input: x, asset price vector, p, momentum vector
% Output: second derivative of Hamilton with respect to p, momentum vectors

    global oneOverGamma;
    val = oneOverGamma * instCov(x);
end
