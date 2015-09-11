function y = hamDp(x, p, mu, sig, gamma)
    y = 1/gamma * (sig*sig*x*x*p + mu*x);
end

