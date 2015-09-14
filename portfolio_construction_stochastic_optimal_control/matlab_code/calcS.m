function [S] = calcS(x, t, T, timeStep, tol, maxIter)
% Input: x, initial asset price vector, t: starting time, T: terminal time
% Output: S: action from WKB approximation
% Note: make sure (T - t) / timeStep is even, i.e. numStep is even

  F = length(x);
  termVal = solveForTermX(x, t, T, timeStep, tol, maxIter)
  zeroVec = zeros(F, 1);

  [xPath, pPath] = generateLfFlow(termVal, zeroVec, t, T, timeStep, tol, maxIter);

  numSteps = size(xPath, 2) - 1;   # number of step equals to number of point minus one.

  # dx <- matrix(nrow = numSteps - 1, ncol = F);
  # dp <- matrix(nrow = numSteps - 1, ncol = F);
  # dpDx <- array(dim = c(numSteps-1, F, F));
  # for(i in 1:numSteps - 1)
  #{
  #dx[i,] <- xPath[i+1,] - xPath[i,];
  #dp[i,] <- pPath[i+1,] - pPath[i,];
  # dpDx[i,,] <- outer(dp[i], dx[i], FUN = "/");
  #}
       
    S0 = 0.0;
    S1 = 0.0;
    for i = 1:(numSteps/2)   # Note the original R code is different here
    
      S0 = S0 + lagr(xPath(:,2*i-1), pPath(:, 2*i-1)) + 4 * lagr(xPath(:, 2*i), pPath(:, 2*i)) + lagr(xPath(:, 2*i+1), pPath(:, 2*i+1));
    # S1 <- S1 + sum(diag(instCov(xPath[2*i-1,]) %*% dpDx[2*i-1,,] + 4 * instCov(xPath[2*i,]) %*% dpDx[2*i,,] + instCov(xPath[2*i+1,]) %*% dpDx[2*i+1,,]));

    end

    S0 = S0 * timeStep / 3;
    S1 = S1 * timeStep / 3;
    
    S = S0 + 0 * S1;
          
end
