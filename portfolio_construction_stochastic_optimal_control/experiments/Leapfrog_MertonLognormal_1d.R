### Leapfrog method for Merton's optimal control of 1 dimensional log normal stochastic process
### 
### H(x,p) = 1/(2 gamma) p C p + 1/gamma p mu x + k/2 mu C^(-1) mu
### Time is from [0,T], divide the interval by N, i.e. h = T/N
### leapfrog iteration is given by
### p_{k-1/2} = p_k + 1/gamma p_k (sig^2 x_k + mu) * h/2
### x_{k-1}   = x_k - 1/gamma (sig^2 x_k^2 p_{k-1/2} + mu x_k) * h (in the implementation, this equation changed)
### p_{k-1}   = p_{k-1/2} + 1/gamma p_{k-1/2}(sig^2 x_{k-1} + mu) * h/2
###




hamDx <- function(x, p, mu, sig, gamma)
{
  hamDx <- 1/gamma * p * (sig*sig*x + mu)
}

hamDp <- function(x, p, mu, sig, gamma)
{
  hamDp <- 1/gamma * (sig*sig*x*x*p + mu*x)
}


xpFlows <- function(period, N, xT, pT, mu, sig, gamma)
{
  flows <- matrix(data=0.0, nrow=N, ncol=2)
  flows[N,1] <- xT
  flows[N,2] <- pT
  
  h <- period/N
  
  for (i in (N-1):1)
  {
    pHalf <- flows[i+1,2] + 0.5 * h * hamDx(flows[i+1,1], flows[i+1,2], mu, sig, gamma)
    #flows[i,1] <- flows[i+1, 1] - h * hamDp(flows[i+1,1], pHalf, mu, sig, gamma)
    flows[i,1] <- (1 - mu*h/2/gamma)/(1 + mu*h/2/gamma) * flows[i+1,1]
    flows[i,2] <- pHalf + 0.5 * h * hamDx(flows[i,1], pHalf, mu, sig, gamma)
  }
  
  return(flows)
}


mu <- 0.05
sig <- 0.19
gamma <- 0.01
xT <- 0.1
pT <- 0
period <- 2
N <- 100

flow <- xpFlows(period, N, xT, pT, mu, sig, gamma)

plot(flow[,1], flow[,2], xlab='position x', ylab='momentum p')

### The exact solution to the 1d lognormal model
t <- seq(from=0, to=period, length.out=(N+1))[-1]
x <- xT * exp(-mu/gamma*(period-t))
p <- 0

### Plot the leapfrog solution and overlapped with the exact solution
plot(1:N, flow[,1], xlab='time t', ylab='position x', col='red')
points(1:N, x, col='blue')















 





