function [val] = lagr(x, p)
  % Lagrange function
  % Input: x, asset price vector, p vector,
  % Output: lagrange value
  % Note: In the R code, it uses term1 - term2, I am not sure if it is what it should be. I am suspecting that the original code ignore the minus sign before According to the paper, this should be term2 - term1, which is what I am using here. this need to be double checked later.

  global corrMatr;
  global oneOverGamma;
  global kappa;
  
  a = driftV(x);
  
  term1 = 0.5 * oneOverGamma * p' * instCov(x) * p;
  term2 = 0.5 * kappa * a' * invInstCov(x) * a;
  val = term2 - term1;

end
