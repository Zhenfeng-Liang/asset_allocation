#_______________________________________________________________________________________________

# Generate the flow of the Hamiltonian system given the terminal condition

#_______________________________________________________________________________________________

generateLfFlow <- function(xT, pT, numFac, hamDx, hamDp, t, T, timeStep, tol, maxIter)
{
  numSteps <- ceiling((T - t) / timeStep);
  h <- 0.5 * T / numSteps;
   
  xFlow <- matrix(data = 0.0, nrow = numFac, ncol = numSteps);
  pFlow <- matrix(data = 0.0, nrow = numFac, ncol = numSteps);
   
  x <- xT;
  p <- pT;
  for(i in numSteps:2)
  {
    xFlow[,i] <- x;
    pFlow[,i] <- p;
    
    err <- 1.0;
    iter <- 0;
    while((err > tol) & (iter <= maxIter))
    {
      pNew <- p + h * hamDx(x, pOld);
      incr <- pNew - pOld;
      err <- sqrt(sum(incr ^ 2));
      pOld <- pNew;
      iter <- iter + 1;
    }
    
    err = 1.0;
    iter = 0;
    while((err > tol) & (iter <= maxIter))
    {
      xNew <- x - h * (hamDp(x, pNew) + hamDp(xOld, pNew));
      incr <- xNew - xOld;
      err <- sqrt(sum(incr ^ 2));
      iter <- iter + 1;
    }
    
    x = xNew;
    p = pNew + h * hamDx(xNew, pNew);
    
    x <- z[1:F];
    p <- z[F:(2*F)];
  }
  xFlow[,1] = x;
  pFlow[,1] = p;
    
  res <- rbind(xFlow, pFlow);
  return(res);
}

#_______________________________________________________________________________________________

oneTimeStep <- function(x, p, hamDx, hamDp, h, tol, maxIter)
{  
  pOld <- p;
  pNew <- p;
  xOld <- x;
  xNew <- x;
  
  err <- 1.0;
  iter <- 0;
  while((err > tol) & (iter <= maxIter))
  {
    pNew <- p + h * hamDx(x, pOld);
    incr <- pNew - pOld;
    err <- sqrt(sum(incr ^ 2));
    pOld <- pNew;
    iter <- iter + 1;
  }
   
  err = 1.0;
  iter = 0;
  while((err > tol) & (iter <= maxIter))
  {
    xNew <- x - h * (hamDp(x, pNew) + hamDp(xOld, pNew));
    incr <- xNew - xOld;
    err <- sqrt(sum(incr ^ 2));
    iter <- iter + 1;
  }
    
  x = xNew;
  p = pNew + h * hamDx(xNew, pNew);
  
  res <- rbind(x, p);
 
  print(res)
  
  return(res);
}