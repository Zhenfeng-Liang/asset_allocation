function y = hamDxp(x, p, rho, mu, sig, gamma)
    y = 1/gamma * (cDx(x, rho, sig)*p + aDx(x, mu));
end