function [val] = ham(x, p)
  % Hamilton function
  % Input: x, asset price vector, p, momentum vector
  % Output: hamiltonian value

  global oneOverGamma;
  global kappa;

  a = driftV(x);
  
  term1 = 0.5 * oneOverGamma * p' * instCov(x) * p;
  term2 = oneOverGamma * p' * a;
  term3 = 0.5 * kappa * a' * invInstCov(x) * a;
  val = term1 + term2 + term3;
  
end
