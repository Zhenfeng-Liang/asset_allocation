function x = newtonRaphsonX(tol, hamilton, xPrev, pPrev, h)
    % Initial value for the iteration
    xOld = xPrev;
    [nrow,~] = size(xPrev);
    % Initial error
    error = 1;
    % Loop through to find implicit solution for x
    while error>tol
        f = xOld + 0.5*h*hamilton.hamDp(xOld, pPrev)'...
            - xPrev + 0.5*h*hamilton.hamDp(xPrev, pPrev)';
        Df = eye(nrow) + 0.5*h*hamilton.hamDxp(xOld, pPrev);
        xNew = double(xOld - inv(Df)*f);
        error = abs(xNew - xOld);
        % symbolic toolbox use rational numbers. Matlab treat it as string
        xOld = double(xNew);
    end
    x = xOld;
end