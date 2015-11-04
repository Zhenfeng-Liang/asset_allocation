function paper1()

    close all;
    clc;
    clear;

    runLN();

end




function runLN()
% run all the lognormal stuff

    display('Start checking lognormal model');
        
    modelParam.modelType = 'LogNormal';
    modelParam.mu = [0.8];
    modelParam.vol = [0.33];
    
    pT = [0];
    
    corrMatr = eye(1);
    gamma = 10.0;
    xCurr = [1.0];
    tCurr = 0;
    T = 10.0;
    timeStep = 0.01;
    tol = 1e-6;

    w0 = 100000;
    numCores = 1;
    utilityType = 'HARA';
    
    turnedOnConsumption = false;
    
    display('Start checking lognormal leapfrog')

    
    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);
    
    utiCalc = UtilityCalculator(gamma, utilityType);
    hamSys = HamiltonianSystem(portCalc, utiCalc);   
    
    wkbSolver = WKBHierarchySolver(hamSys, numCores);
    
    xT = wkbSolver.solveForTermX(xCurr, tCurr, T, timeStep, tol);
    
    [xFlow, pFlow] = wkbSolver.generateLfFlow(xT, pT, tCurr, T, timeStep, tol);
    
    t = tCurr:timeStep:T;    
    x_exact_1 = xT(1) * exp(-modelParam.mu(1) / gamma * (T - t));
 
    display('Plotting lognomal exact xFlow and leapfrog flow on the same graph. They should be on top of each other') 
    
    figure;
    
    subplot(1,2,1);
    plot(t, xFlow(1,:), 'Color', 'red');
    hold on;
    plot(t, x_exact_1, 'Color', 'blue');
    
    title('Asset 1: Lognormal wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label    
    
    axes('Position',[.17 .6 .15 .15])
    box on
    
    diff = xFlow(1,:) - x_exact_1;
    plot(t, diff)
    
    title('Difference');
    xlabel('time') % x-axis label
    ylabel('difference') % y-axis label    

    
    %figure;
    subplot(1,2,2);
    plot(t, pFlow(1,:), 'Color', 'red');
    hold on;
    plot(t, 0, 'Color', 'blue');
    
    title(['Asset 2: Lognormal wkb approximate pFlow and exact ' ...
           'pFlow, should be zero'])
    xlabel('time') % x-axis label
    ylabel('pFlow') % y-axis label
    
    
    display(['Calculating wkb strategy approximation and exact ' ...
             'strategy'])
    
    
    wkbStrategy = wkbOptimizer(modelParam, corrMatr, gamma, xCurr, tCurr, ...
                        T, timeStep, tol, w0, utilityType, ...
                        turnedOnConsumption, numCores)

    exactStrategy = 1.0 / utiCalc.Au(w0) * portCalc.invInstCov(xCurr) * model.driftV(xCurr)
    
    diffNorm = norm(wkbStrategy - exactStrategy)

    
    display('Finished checking lognormal model')

end

