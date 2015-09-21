function [phi] = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter)
% Input: xCurr: current asset price vector, tCurr: current time, T: terminal time, timeStep: time interval between points, tol: tolerance, maxIter: maximum iteration number
% Output: phi: asset allocation weight vector
% Note: make sure (T - tCurr) / timeStep is even, i.e. numStep is odd.

    phi = invInstCov(xCurr) * driftV(xCurr) + calcNablaS(xCurr, tCurr, T, timeStep, tol, maxIter);
end
    

