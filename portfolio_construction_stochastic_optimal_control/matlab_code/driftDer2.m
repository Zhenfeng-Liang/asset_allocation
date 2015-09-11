function [der2] = driftDer2(x, i, j)
% Input: x, asset price vector, with respect to ith and jth asset
% Output: second derivative of drift vector with respect to ith and jth assets

  # mean reverting and lognormal models
  F = length(x);
  der2 = zeros(1, F);    # check the dimension here later. I am suspecting this shoule be F x 1
  
end

