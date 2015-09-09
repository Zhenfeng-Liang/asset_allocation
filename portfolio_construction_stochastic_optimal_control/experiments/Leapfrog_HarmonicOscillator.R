### leapfrog harmonic oscillator
###
### Hamiltonian is H = 0.5 * p^2 + 0.5 * x^2
### Time is from [0, T], divide the interval by N, i.e. h = T/N
### leapfrog interation is given by
### p_{k+1/2} = p_k - h/2 * x_k
### x_{k+1}   = x_k + h * p_{k+1/2}
### p_{k+1}   = p_{k+1/2} - h/2 * x_{k+1}
###

hamDx <- function(x)
{
  return(x)
}

hamDp <- function(p)
{
  return(p)
}

xpFlows <- function(Period, N)
{
  flows <- matrix(data=0.0, nrow=N, ncol=2)
  flows[1,1] = 0
  flows[1,2] = 1
  h <- Period/N
  
  for (i in 2:N)
  {
    pHalf <- flows[i-1, 2] - 0.5*h*hamDx(flows[i-1,1])
    flows[i,1] <- flows[i-1,1] + h*hamDp(pHalf)
    flows[i,2] <- pHalf - 0.5*h*hamDx(flows[i,1])
  }
  
  return(flows)
}

period <- 7
N <- 200
flow <- xpFlows(period, N)
plot(flow[,1],flow[,2],xlab='position x', ylab='momentum p', main='Leapfrog method for simple harmonic oscilator')

