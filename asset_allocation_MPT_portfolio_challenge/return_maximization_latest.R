
variance <- function(w)
{
  # daily variance
  tmp <- t(w)%*%cov_matrix%*%w
  
  return(as.numeric(tmp))
}

portRet <- function(w)
{
  pRet <- as.numeric(w%*%ret_period)
  
  return (1 / pRet)
}

sum_of_weight <- function(w)
{
  tmp <- sum(w)
  return(tmp)
}

constraint <- function(w)
{
   sum_of_weight <- sum(w)
   var <- variance(w)
   #var <- as.numeric(t(w)%*%cov_matrix%*%w)
   return(c(sum_of_weight, var))
}


############################################ Main script #####################################################
require(quantstrat)
require(Rsolnp)

root_path = Sys.getenv('PROJ_ROOT')
data_dir = paste(root_path,'data_semi_annual',sep='')  

read_data_locally = TRUE

start = "2014-12-01"
end = "2015-05-31"

total_capital <- 1000000

ticker_universe <- c("CAT", "DYN", "VRTX", "UTHR", "HUM", "ACT", "CXO", 
                     "SCTY", "AMGN", "GILD", "CELG", "PFE", "JNJ", "UNH", "LLY", 
                     "ABT", "BDX", "SYK", "ISRG", "MON", "DXCM", "GEVA",         # new stocks
                     "ISIS", "MIC", "CI")                                        # new stocks

stock_num <- length(ticker_universe)

if(read_data_locally) {
  getSymbols(ticker_universe, src='rda',dir=data_dir, extension='RData', from=start, to=end) 
  getSymbols("RUI", src='rda',dir=data_dir, extension='RData', from=start, to=end) 
}else {
  getSymbols(ticker_universe, from = start, to = end)
  getSymbols("^RUI", from = start, to = end)
  saveSymbols(ticker_universe, data_dir)
  saveSymbols("RUI", data_dir)
}

rui_xts <- get("RUI")
rui_ret <- dailyReturn(rui_xts)
rui_daily_var <- as.numeric(var(rui_ret))

ret_matrix <- rui_ret
ret_period <- c()
for(ticker in ticker_universe)
{
  temp_xts <- get(ticker)
  temp_ret <- dailyReturn(temp_xts)
  ret_matrix <- cbind(ret_matrix, temp_ret)
  
  end_price <- as.numeric(Cl(tail(temp_xts, n = 1)))
  start_price <- as.numeric(Cl(temp_xts[1,]))
  ret_period <- c(ret_period, (end_price - start_price)/start_price)
}

# daily
ret_matrix <- ret_matrix[,-1]
cov_matrix <- cov(ret_matrix)

dim(ret_period) <- c(stock_num, 1)

w0 <- rep(1/stock_num, stock_num)

lowerB <- rep(0.02, stock_num)
upperB <- rep(0.08, stock_num)

opWeights <- solnp(w0, fun = portRet, eqfun = sum_of_weight, eqB = 1.0, 
                   ineqfun = variance, ineqUB = 1.5 * rui_daily_var, ineqLB = 1.0 * rui_daily_var, 
                   LB=lowerB, UB=upperB)

#opWeights <- solnp(w0, fun = portRet, eqfun = constraint, eqB = c(1.0, 2.5 * rui_daily_var), 
#                   LB=lowerB, UB=upperB)

rui_ret <- (as.numeric(Cl(tail(rui_xts, n=1))) - as.numeric(Cl(rui_xts[1,]))) / as.numeric(Cl(rui_xts[1,]))

cat("Russell 1000's return is ", rui_ret, " Portfolio return is ", 1 / tail(opWeights$values, n=1))

allocation <- opWeights$pars * total_capital
View(allocation)
