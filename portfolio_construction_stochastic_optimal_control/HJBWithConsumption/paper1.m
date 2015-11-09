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
    modelParam.mu = [0.3];
    modelParam.vol = [0.33];
    
    pT = [0];
    
    corrMatr = eye(1);
    gamma = 10.0;
    xCurr = [1.0];
    tCurr = 0;
    T = 1.0;
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
    legend('Numerical','Exact');
    
    
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
    legend('Numerical','Exact');

    
    btST = 0;           
    btET = 1.0;               
    rebTS = 0.01;


    simulator = ModelEvolver();
    simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                     model);

    F = size(simData, 1);
    m = size(simData, 2);
        
    constr = Constraint(0.3, -0.2, false);
    bte = BtEngine(btST, btET, rebTS, constr);

    [wVec, phiMat, cVec] = bte.runBackTest(simData, wkbSolver, w0, ...
                                           turnedOnConsumption);

    exactStrategy = zeros(F, m-1);
    
    for i = 1:(m-1)        
        exactStrategy(:,i) = 1.0 / utiCalc.Au(wVec(1,i)) * portCalc.invInstCov(simData(:,i)) ...
            * model.driftV(simData(:,i));        
    end

    t = btST:rebTS:(btET-rebTS);    

    figure;
    
    plot(t, phiMat(1,:), 'Color', 'red');
    hold on;
    plot(t, exactStrategy(1,:), 'Color', 'blue');
    
    title('LogNormal numerical and exact strategies');
    xlabel('rebalance time') % x-axis label
    ylabel('phi*') % y-axis label    
    legend('Numerical','Exact');
    
    axes('Position',[.17 .2 .2 .2])
    box on
    
    diff = phiMat(1,:) - exactStrategy;
    plot(t, diff)
    
    title('Difference');
    xlabel('rebalance time') % x-axis label
    ylabel('difference') % y-axis label    

    
    t = btST:rebTS:btET; 
    pos = [phiMat(1,:), 0];
    figure;
    [ax, h1, h2] = plotyy(t, simData(1,:), t, pos, 'plot', 'stairs');    
    set(h1, 'linewidth', 2, 'color', 'red');
    set(h2,  'linewidth', 2, 'color', 'blue');
    
    title('Simulated Stock Price and WKB approximate position', ...
          'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20) % x-axis label
    ylabel(ax(1), 'simulated stock price($)', 'color', 'red', 'FontSize', 20);
    ylabel(ax(2), 'position', 'color', 'blue', 'FontSize', 20);
    
    h = legend('simulated stock price','wkb strategy', 'Location', 'southeast');
    set(h, 'FontSize', 16);

    % bar(t, [0,posChg,0]);
    
    
    figure;
    plot(t, wVec(1,:), 'Color', 'red');    
    title('Cumulative PnL');
    xlabel('rebalance time') % x-axis label
    ylabel('PnL') % y-axis label    
    
end

