function y = hamDp(x, p, rho, mu, sig, gamma)
    y = 1/gamma * (c(x, rho, sig)*p + a(x, mu));
end

