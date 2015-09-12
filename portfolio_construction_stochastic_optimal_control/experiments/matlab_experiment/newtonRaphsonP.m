function p = newtonRaphsonP(tol, hamDx_h, hamDxp_h, xPrev, pPrev, h, rho, mu, sig, gamma)
    % Initial value for the iteration
    pOld = pPrev;
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = pOld - pPrev - 0.5*h*feval(hamDx_h, xPrev, pOld, rho, mu, sig, gamma);
        Df = 1 - 0.5*h*feval(hamDxp_h, xPrev, pOld, rho, mu, sig, gamma);
        pNew = pOld - f/Df;
        error = abs(pNew - pOld);
        pOld = pNew;
    end
    p = pOld;
end