function y = c(x, rho, sig)
    %y = rho*sig*sig*((0.05-x)^2);
    y = rho*sig*sig*x*x;
end