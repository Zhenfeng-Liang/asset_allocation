#_______________________________________________________________________________________________

# Model definition

#_______________________________________________________________________________________________

driftV <- function(x)
{
  # mean reverting model
  mu <- c(0.4, -1.3, -2.5, 0.1, 5.9);
  lambda <- c(21.0, 13.2, 7.4, 24.1, 18.1);
  d <- lambda * (mu - x);
  
  # lognormal model
  #mu <- c(0.4, 0.12, 0.05, 0.16, 0.02);
  #d <- mu * x;                     
  
  return(d);
  
}

#_______________________________________________________________________________________________

diffV <- function(x)
{  
  # mean reverting model
  nVols <- c(13.2, 23.0, 21.1, 7.5, 18.7);
  d <- diag(nVols);
  
  # lognormal model
  #lnVols <- c(0.33, 0.26, 0.19, 0.22, 0.31);
  #d <- diag(lnVols * x);    
  
  return(d);
}

#_______________________________________________________________________________________________

driftDer <- function(x, i)
{ 
  # mean reverting model
  lambda <- c(21.0, 13.2, 7.4, 24.1, 18.1);
  F <- length(lambda);
  der <- vector(mode = "double", F);
  der[i] <- -lambda[i];
  
  # lognormal model
  #mu <- c(0.4, 0.12, 0.05, 0.16, 0.02);
  #F <- length(mu);
  #der <- vector(mode = "double", F);
  #der[i] <- mu[i];
  
  return(der);
}

#_______________________________________________________________________________________________

diffDer <- function(x, i)
{ 
  # mean reverting model
  F <- 5;
  der <- matrix(data = 0.0, nrow = F, ncol = F);
  
  # lognormal model
  #lnVols <- c(0.33, 0.26, 0.19, 0.22, 0.31);
  #F <- length(lnVols);
  #vecD <- vector(mode = "double", F);
  #vecD[i] <- lnVols[i];
  #der <- diag(vecD);
  
  return(der);
  
}

#_______________________________________________________________________________________________

driftDer2 <- function(x, i, j)
{ 
  # mean reverting and lognormal models
  F <- 5;
  der <- vector(mode = "double", F);   
  
  return(der);
}

#_______________________________________________________________________________________________

diffDer2 <- function(x, i, j)
{  
  # mean reverting and lognormal models
  F <- 5;
  der <- matrix(data = 0.0, nrow = F, ncol = F);   
  
  return(der);
}