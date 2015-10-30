function checkLogNormal()
% This function is for test purpose. For lognomal model, it plots the wkb
% approximate xFlow and the exact xFlow. They should be on the top of each
% other. It also calculates the norm of difference between wkb approximate
% strategy and the exact strategy, which should be very small.   

    display('Start checking lognormal model');
    
    
    modelParam.modelType = 'LogNormal';
    modelParam.mu = [0.8; 1.3; 0.9];
    modelParam.vol = [0.33; 0.16; 0.45];
    
    corrMatr = eye(3);
    gamma = 10.0;
    xCurr = [1.0; 0.8; 1.5];
    tCurr = 0;
    T = 0.25;
    timeStep = 0.01;
    tol = 0.0001;

    w0 = 100000;
    numCores = 2;
    utilityType = 'CRRA';
    
    turnedOnConsumption = true;
    
    display('Start checking lognormal leapfrog')

    
    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);
    
    utiCalc = UtilityCalculator(gamma, utilityType);
    hamSys = HamiltonianSystem(portCalc, utiCalc);   
    
    wkbSolver = WKBHierarchySolver(hamSys, numCores);
    
    xT = wkbSolver.solveForTermX(xCurr, tCurr, T, timeStep, tol);
    
    [xFlow, pFlow] = wkbSolver.generateLfFlow(xT, [0;0;0], tCurr, T, timeStep, tol);
    
    t = tCurr:timeStep:T;    
    x_exact_1 = xT(1) * exp(-modelParam.mu(1) / gamma * (T - t));
    x_exact_2 = xT(2) * exp(-modelParam.mu(2) / gamma * (T - t));
    x_exact_3 = xT(3) * exp(-modelParam.mu(3) / gamma * (T - t));
 
    display('Plotting lognomal exact xFlow and leapfrog flow on the same graph. They should be on top of each other') 
    
    figure;
    
    subplot(2,2,1);
    plot(t, xFlow(1,:), 'Color', 'red');
    hold on;
    plot(t, x_exact_1, 'Color', 'blue');
    
    title('Asset 1: Lognormal wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label    
    
    subplot(2,2,2);
    plot(t, xFlow(2,:), 'Color', 'red');
    hold on;
    plot(t, x_exact_2, 'Color', 'blue');
    
    title('Asset 2: Lognormal wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label
    
    subplot(2,2,3);
    plot(t, xFlow(3,:), 'Color', 'red');
    hold on;
    plot(t, x_exact_3, 'Color', 'blue');

    title('Asset 3: Lognormal wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label
    
    display(['Calculating wkb strategy approximation and exact ' ...
             'strategy'])
    
    wkbStrategy = wkbOptimizer(modelParam, corrMatr, gamma, xCurr, tCurr, ...
                        T, timeStep, tol, w0, utilityType, ...
                        turnedOnConsumption, numCores)

    exactStrategy = 1.0 / utiCalc.Au(w0) * portCalc.invInstCov(xCurr) * model.driftV(xCurr)
    
    diffNorm = norm(wkbStrategy - exactStrategy)

    
    display('Finished checking lognormal model')

end

