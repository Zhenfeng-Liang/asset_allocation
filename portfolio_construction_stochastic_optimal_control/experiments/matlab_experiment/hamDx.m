function y = hamDx(x, p, mu, sig, gamma) 
    y = 1/gamma * p * (sig*sig*x + mu);
end