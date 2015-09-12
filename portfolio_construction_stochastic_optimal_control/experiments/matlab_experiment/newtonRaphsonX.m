function x = newtonRaphsonX(tol, hamDp_h, hamDxp_h, xPrev, pPrev, h, rho, mu, sig, gamma)
    % Initial value for the iteration
    xOld = xPrev;
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = xOld + 0.5*h*feval(hamDp_h, xOld, pPrev, rho, mu, sig, gamma)...
            - xPrev + 0.5*h*feval(hamDp_h, xPrev, pPrev, rho, mu, sig, gamma);
        Df = 1 + 0.5*h*feval(hamDxp_h, xOld, pPrev, rho, mu, sig, gamma);
        xNew = xOld - f/Df;
        error = abs(xNew - xOld);
        xOld = xNew;
    end
    x = xOld;
end