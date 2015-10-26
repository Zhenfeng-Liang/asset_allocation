function [strategy] = wkbOptimizer(modelParam, corrMatr, gamma, xCurr, tCurr, T, timeStep, tol)
% Assumed exponential utility
% Input:
%       1: modelParam, struct, should include:
%               modelParam.modelType: string, 'MeanReverting or LogNormal')
%               modelParam.mu: mean column vector in respective models
%               modelParam.vol: diffusion coeff column vector for every Brownian Motion
%               modelParam.lambda: optional, for MR model, column vector,
%                                  mean reverting speed for each asset
%
%       2: corrMatr: correlation matrix between brownian motion. Dim: p*p
%
%       3: gamma: parameters for exponential utility function
%
%       4: xCurr: column vector, current asset price
%
%       5: T: terminal time
%
%       6: timeStep: time step for leapfrog algorithm
%
%       7: tol: tolerance for Newton methods

    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);
    hamSys = HamiltonianSystem(portCalc, gamma);
    wkbSolver = WKBHierarchySolver(hamSys);
    
    display(['Start optimizing portfolio under ', modelParam.modelType, ...
            ' model']);
    tic
    %strategy = 1.0 / gamma * wkbSolver.optimalControlStrategy(xCurr, ...
    %                                                  tCurr, T, timeStep, tol);
    strategy = wkbSolver.optimalControlStrategy(xCurr, tCurr, T, timeStep, tol);
    toc
    
    display(['Finished optimizing portfolio under ', modelParam.modelType, ...
            ' model']);

end

