function p = newtonRaphsonP(tol, hamilton, xPrev, pPrev, h)
    % Initial value for the iteration
    pOld = pPrev;
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = pOld - pPrev - 0.5*h*hamilton.hamDx(xPrev, pOld);
        Df = 1 - 0.5*h*hamilton.hamDxp(xPrev, pOld);
        pNew = pOld - f/Df;
        error = abs(pNew - pOld);
        pOld = pNew;
    end
    p = pOld;
end