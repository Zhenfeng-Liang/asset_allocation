# The solver is to minimize the obj func, so, in order to maximize the function, we use
# (total_capital - portfolio_position) / total_capital to calculate obj value
fn1 <- function(x){

  portfolio_position <- 0;
  for(i in 1:stock_num)
  {
    tmpxts <- get(ticker_universe[i])
    start_closed_price <- as.numeric(Cl(tmpxts[1,]))
    len <- length(index(tmpxts))
    end_closed_price <- as.numeric(Cl(tmpxts[len,]))

    initial_quantity <- x[i] / start_closed_price
    
    portfolio_position <- portfolio_position + initial_quantity * end_closed_price
  }
  
  ret <- total_capital / (portfolio_position - total_capital) 
  #ret <- (total_capital - portfolio_position) / total_capital
  
  return (ret)
}

# Two constraints:
# 1: sum of x equals to total_capital
# 2: standard deviation equals to RUI standard equation
eqn1 <- function(x){
  
  initial_position <- sum(x)
  
  initial_quantity <- c()
  for(i in 1:stock_num) 
  {
    tmpxts <- get(ticker_universe[i])
    start_closed_price <- as.numeric(Cl(tmpxts[1,]))
    initial_quantity <- c(initial_quantity, x[i] / start_closed_price)
  }
  
  len <- length(index(get("RUI")))
  portfolio_position <- c()
  for(j in 1:len)
  {
    current_portfolio_position <- 0
    for(i in 1:stock_num)
    {
      tmpxts <- get(ticker_universe[i])
      current_stock_price <- as.numeric(Cl(tmpxts[j,]))
      current_portfolio_position <- current_portfolio_position + current_stock_price * initial_quantity[i]
    }
    portfolio_position <- c(portfolio_position, current_portfolio_position)
  }
  
  sd_port <- sd(portfolio_position)
  
  return(c(initial_position, sd_port))
}


############################################ Main script #####################################################
require(quantstrat)
require(Rsolnp)

root_path = Sys.getenv('PROJ_ROOT')
data_dir = paste(root_path,'data_two_month',sep='')  

read_data_locally = TRUE

start = "2015-04-01"
end = "2015-05-31"

total_capital <- 1000000

ticker_universe <- c("CAT", "DYN", "VRTX", "UTHR", "HUM", "ACT", "CXO", 
                    "SCTY", "AMGN", "GILD", "CELG", "PFE", "JNJ", "UNH", "LLY")

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
rui_cl_price <- as.numeric(Cl(rui_xts))
rui_initial_quantity <- total_capital / rui_cl_price[1]
rui_position <- sapply(rui_cl_price, function(x) x * rui_initial_quantity)
rui_sd <- sd(rui_position)

rui_ret <- (tail(rui_cl_price, n=1) - rui_cl_price[1]) / rui_cl_price[1]

initial_each_stock_allocation <- total_capital / stock_num
initial_guess <- rep(50000, stock_num)

lowerB <- rep(0, stock_num)
upperB <- rep(total_capital, stock_num)

result=solnp(initial_guess, fun = fn1, eqfun = eqn1, eqB = c(total_capital, rui_sd), 
             LB=lowerB, UB=upperB)

cat("Russell 1000's return is ", rui_ret, " Portfolio return is ", 1 / tail(result$values, n=1))

print(result$elapsed)

allocated_position <- result$pars
View(allocated_position)
