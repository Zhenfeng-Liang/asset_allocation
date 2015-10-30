function [strategy] = wkbOptimizer(modelParam, corrMatr, gamma, ...
                                   xCurr, tCurr, T, timeStep, tol, ...
                                   w0, utilityType, turnedOnConsumption, ...
                                   numCores)
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
%
%       8: w0: initial wealth
%
%       9: utilityType: utility function type. Right now, only CRRA
%       works for consumption part.
%
%       10: turnedOnConsumption: boolean, true for consumption on,
%       otherwise off
%
%       11: numCores: number of cores to run the program, only
%       works for consumption part.
    
    % open a log file
    fid=fopen('./log.txt','wt');
    
    % optimizer set up
    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);
    
    utiCalc = UtilityCalculator(gamma, utilityType);
    
    hamSys = HamiltonianSystem(portCalc, utiCalc);
    
    wkbSolver = WKBHierarchySolver(hamSys, numCores);
    
    display(['Start optimizing portfolio under ', modelParam.modelType, ...
            ' model']);
    tic;
    
    % standard out to show optimizer set up
    numCores
    xCurr
    tCurr
    T
    utilityType
    
    % Run the strategy
    [strategy, consumingStrategy] = wkbSolver.optimalControlStrategy(xCurr, tCurr, T, ...
                                                      timeStep, tol, w0, turnedOnConsumption)
    toc
    ttime = toc;
    
    display(['Finished optimizing portfolio under ', modelParam.modelType, ...
            ' model']);
    
    % Writing the result into log file
    fprintf(fid, 'consumption turned on: %f\n', turnedOnConsumption);
    fprintf(fid,'Number of processors (labs) used was: %d\n', numCores);
    fprintf(fid,'Strategy is\n');
    fprintf(fid, [repmat('%f\t', 1, size(strategy, 2)) '\n'], strategy');
    fprintf(fid,'Time to complete the computation was: %6.6f\n', ttime);
    
    fclose(fid);
end

