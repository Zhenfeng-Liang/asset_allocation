function [d] = driftV(x)
% Given x, asset price, output asset dynamics drift
% Right now, I am hard coding two dimensions drift and lambda
% Return: drift vector

  # mean reverting model
  mu = [0.4, -1.3];
  lambda = [21.0, 13.2];
  d = lambda .* (mu - x);
  
end
