library(MASS)
library(expm)

#_______________________________________________________________________________________________

# F-factor optimal portfolio problem

#_______________________________________________________________________________________________

wkbOptPortf <- function(xCurr, tCurr, T, driftV, diffV, driftDer, diffDer, driftDer2, diffDer2, corrMatr, gamma, timeStep, tol, maxIter)
{
  F <- ncol(corrMatr);
  oneOverGamma <- 1 /gamma;
  kappa <- oneOverGamma - 1;
  
  # Auxilairy functions
  instCov <- function(x)
  {    
    dofx <- diffV(x);
    covx <- t(dofx) %*% corrMatr %*% dofx;
    
    return(covx);
  }
  
  invInstCov <- function(x)
  {    
    invCovx <- solve(instCov(x));
    
    return(invCovx);
  }
  
  instCovDer <- function(x, i)
  { 
    dofx <- diffV(x);
    ddofx <- diffDer(x, i);
    derx <- t(ddofx) %*% corrMatr %*% dofx + t(dofx) %*% corrMatr %*% ddofx;
      
    return(derx);
  }
  
  invInstCovDer <- function(x, i)
  {    
    iCovx <- invInstCov(x);
    derx <- -iCovx %*% instCovDer(x, i) %*% iCovx;
    
    return(derx);
  }
  
  instCovDer2 <- function(x, i, j)
  { 
    dofx <- diffV(x);
    ddofxi <- diffDer(x, i);
    ddofxj <- diffDer(x, j);
    derx <- t(diffDer2(x, i, j)) %*% corrMatr %*% dofx + t(ddofxi) %*% corrMatr %*% ddofxj 
                                + t(ddofxj) %*% corrMatr %*% ddofxi +t(dofx) %*% corrMatr %*% diffDer2(x, i, j);
    
    return(derx);
  }
  
  invInstCovDer2 <- function(x, i, j)
  {    
    iCovx <- invInstCov(x);
    iCovxdi <- instCovDer(x, i);
    iCovxdj <- instCovDer(x, j)
    derx <- iCovx %*% (iCovxdi %*% iCovx %*% iCovxdj + iCovxdj %*% iCovx %*% iCovxdi - instCovDer2(x, i, j)) %*% iCovx;
    
    return(derx);
  }
  
  # Lagrange function
  lagr <- function(x, p)
  {
    d <- driftV(x);
      
    term1 <- 0.5 * oneOverGamma * t(p) %*% instCov(x) %*% p;
    term2 <- 0.5 * kappa * t(d) %*% invInstCov(x) %*% d;
    val <- term1 - term2;
     
    return(val);
  }
  
  # Hamilton function
  ham <- function(x, p)
  {
    d <- driftV(x);
    
    term1 <- 0.5 * oneOverGamma * t(p) %*% instCov(x) %*% p;
    term2 <- oneOverGamma * t(p) %*% d;
    term3 <- 0.5 * kappa * t(d) %*% invInstCov(x) %*% d;
    val <- term1 + term2 + term3;
    
    return(val);
  }
  
  hamDx <- function(x, p)
  {  
    d <- driftV(x);  
    
    term1 <- vector(mode = "double", length = F);
    term2 <- vector(mode = "double", length = F);
    term3 <- vector(mode = "double", length = F);
    
    for(i in 1:F)
    {
      term1[i] <- 0.5 * oneOverGamma * t(p) %*% instCovDer(x, i) %*% p;
      term2[i] <- oneOverGamma * t(p) %*% driftDer(x, i);
      term3[i] <- kappa * t(d) %*% (invInstCov(x) %*% driftDer(x, i) + 0.5 * invInstCovDer(x, i) %*% d);
    }
    
    val <- term1 + term2 + term3;
    
    return(val);
  }
  
  hamDp <- function(x, p)
  {
    val <- oneOverGamma * (instCov(x) %*% p + driftV(x));
          
    return(val);
  }
  
  hamDxx <- function(x, p)
  {    
    d <- driftV(x);  
    
    term1 <- matrix(nrow = F, ncol = F);
    term2 <- matrix(nrow = F, ncol = F);
    term3 <- matrix(nrow = F, ncol = F);
    
    for(i in 1:F)
    {
      for(j in 1:F)    
      {
        term1[i,j] <- 0.5 * oneOverGamma * t(p) %*% instCovDer2(x, i, j) %*% p;
        term2[i,j] <- oneOverGamma * t(p) %*% driftDer2(x, i, j);
        term3[i,j] <- kappa * (t(d) %*% invInstCov(x) %*% driftDer2(x, i, j) 
                               + t(driftDer(x, i)) %*% invInstCov(x) %*% driftDer(x, j)
                               + t(d) %*% invInstCovDer(x, i) %*% driftDer(x, j)
                               + t(d) %*% invInstCovDer(x, j) %*% driftDer(x, i)
                               + 0.5 * t(d) %*% invInstCovDer2(x, i, j) %*% d);
      }
    }
    
    val <- term1 + term2 + term3;
    
    return(val);
  }
  
  hamDxp <- function(x, p)
  {  
    val <- matrix(nrow = F, ncol = F);
    
    for(i in 1:F)
    {
      val[i,] <- oneOverGamma * (instCovDer(x, i) %*% p + driftDer(x, i));
    }
    
    return(val);
  }
  
  hamDpp <- function(x, p)
  {  
    val <- oneOverGamma * instCov(x);
    
    return(val);
  }
   
  # Solve the terminal value problem for Hamilton's equations
  generateLfFlowOptPortf <- function(xT, pT, t)
  {     
    flow <- generateLfFlow(xT, pT, F, hamDx, hamDp, t, T, timeStep, tol, maxIter);
       
    return(flow);
  }
  
  # Invert the Hamiltonian flow (zero terminal momentum)  
  solveForTermX <- function(x, t)
  {
    # setup the system of approximate variational equations
    zeroVec <- rep(0.0, F);
    termVal <- rbind(diag(rep(1.0, F)), diag(rep(0.0, F)));
    
    xT <- x;    
    err <- 1.0;
    ctr <- 1;    
    while(err >= tol && ctr <= maxIter)
    {
      Flow <- generateLfFlowOptPortf(xT, zeroVec, t);
      Q <- hamDxp(xT, zeroVec);
      R <- hamDpp(xT, zeroVec);
      U <- hamDxx(xT, zeroVec);
      apprVarSol <- expm((t-T) * rbind(cbind(Q, R), cbind(-U, -Q))) %*% termVal;
      invDF <- solve(apprVarSol[1:F,]);
      z <- xT - invDF %*% (Flow[1,1:F] - x);
      err <- norm(z - xT);
      xT <- z;
      ctr = ctr + 1;
    }
    print("DONE INVERTING THE FLOW")
    print(xT)
    
    return(xT);    
  }
  
  calcS <- function(x, t)
  {
    termVal <- solveForTermX(x, t);
    zeroVec <- rep(0.0, F);
    flow <- generateLfFlowOptPortf(termVal, zeroVec, t);
    numSteps <- nrow(flow);
    xPath <- flow[,1:F];
    pPath <- flow[,(F+1):(2*F)];
    dx <- matrix(nrow = numSteps - 1, ncol = F);
    dp <- matrix(nrow = numSteps - 1, ncol = F);
    dpDx <- array(dim = c(numSteps-1, F, F));
    for(i in 1:numSteps - 1)
    {
      dx[i,] <- xPath[i+1,] - xPath[i,];
      dp[i,] <- pPath[i+1,] - pPath[i,];
      dpDx[i,,] <- outer(dp[i], dx[i], FUN = "/");
    }
       
    S0 <- 0.0;
    S1 <- 0.0;
    for(i in 1:(numSteps/2 - 1))
    {
      S0 <- S0 + lagr(xPath[2*i-1,], pPath[2*i-1,]) + 4 * lagr(xPath[2*i,], pPath[2*i,]) + lagr(xPath[2*i+1,], pPath[2*i+1,]);
      S1 <- S1 + sum(diag(instCov(xPath[2*i-1,]) %*% dpDx[2*i-1,,] + 4 * instCov(xPath[2*i,]) %*% dpDx[2*i,,] + instCov(xPath[2*i+1,]) %*% dpDx[2*i+1,,]));
    }
    S0 = S0 * timeStep / 3;
    S1 = S1 * timeStep / 3;
    
    S <- S0 + 0 * S1;
          
    return(S);
  }
   
  calcNablaS <- function(x, t)
  {
    S <- calcS(x, t);
    nablaS <- vector(mode = "double", length = F);
    eps <- 1e-5;
    for(i in 1:F)
    {
      xEps <- x;
      xEps[i] <- xEps[i] + eps;
      nablaS[i] = (calcS(xEps, t) - S) / eps;
    }
    
    return(nablaS);
  }
 
  
  # Optimal control (without the utility dependent part)
  #optimalControl <- calcNablaS(xCurr, tCurr);
  #optimalControl <- invInstCov(xCurr) %*% driftV(xCurr);
  optimalControl <- invInstCov(xCurr) %*% driftV(xCurr) + calcNablaS(xCurr, tCurr);
  
  return(optimalControl);
}