function x = newtonRaphsonX(tol, hamilton, xPrev, pPrev, h)
    % Initial value for the iteration
    xOld = xPrev;
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = xOld + 0.5*h*hamilton.hamDp(xOld, pPrev)...
            - xPrev + 0.5*h*hamilton.hamDp(xPrev, pPrev);
        Df = 1 + 0.5*h*hamilton.hamDxp(xOld, pPrev);
        xNew = xOld - f/Df;
        error = abs(xNew - xOld);
        xOld = xNew;
    end
    x = xOld;
end