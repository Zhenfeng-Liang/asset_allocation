function p = newtonRaphsonP(tol, hamilton, xPrev, pPrev, h)
    % Initial value for the iteration
    pOld = pPrev;
    [nrow,~] = size(xPrev);
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = pOld - pPrev - 0.5*h*(hamilton.hamDx(xPrev, pOld))';
        Df = eye(nrow) - 0.5*h*hamilton.hamDxp(xPrev, pOld);
        pNew = double(pOld - inv(Df)*f);
        error = norm(pNew - pOld);
        pOld = double(pNew);
    end
    p = pOld;
end