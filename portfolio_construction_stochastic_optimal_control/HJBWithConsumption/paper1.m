function paper1()

    outdir = '11.19.LR';
    %runLN(outdir);
    %runMR(outdir);
    %runCIR(outdir);

    numPaths = 100;
    
    modelParam.modelType = 'MeanReverting';
    modelParam.mu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5];
    modelParam.vol = [0.1; 0.16; 0.3; 0.52; 0.14; 0.5; 1.0; 0.3; 0.5; 0.8];
    modelParam.lambda = [2.10; 1.32; 1.10; 1.24; 1.56; 0.6; 1.9; 2.3; 1.05; 0.8];
    xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1];

    runMonteCarlo(numPaths, modelParam, xCurr);


    modelParam2.modelType = 'CIR';
    modelParam2.mu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5];
    modelParam2.vol = [0.1; 0.16; 0.3; 0.52; 0.14; 0.5; 1.0; 0.3; 0.5; 0.8];
    modelParam2.vol = modelParam.vol ./ sqrt(modelParam.mu); 
    modelParam2.lambda = [2.10; 1.32; 1.10; 1.24; 1.56; 0.6; 1.9; 2.3; 1.05; 0.8];
    xCurr2 = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1];

    runMonteCarlo(numPaths, modelParam2, xCurr2);

end


function runLN(outdir)
% run all the lognormal stuff
    
    tic
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
    plot(t, xFlow(1,:), 'Color', 'red', 'linewidth', 2);
    hold on;
    plot(t, x_exact_1, 'Color', 'blue', 'linewidth', 2);
    
    title('LN numerical xFlow and exact xFlow', ...
          'FontSize', 20);
    xlabel('time', 'FontSize', 20) % x-axis label
    ylabel('xFlow', 'FontSize', 20) % y-axis label    
    h = legend('Numerical','Exact');
    set(h, 'FontSize', 16);
    
    axes('Position',[.17 .6 .15 .15])
    box on
    
    diff = xFlow(1,:) - x_exact_1;
    plot(t, diff)
    
    title('Difference', 'FontSize', 20);
    xlabel('time', 'FontSize', 20) % x-axis label
    ylabel('difference','FontSize', 20) % y-axis label    

    
    %figure;
    subplot(1,2,2);
    plot(t, pFlow(1,:), 'Color', 'red', 'linewidth', 2);
    hold on;
    plot(t, t * 0, 'Color', 'blue', 'linewidth', 2);
    
    title(['LN numerical pFlow and exact ' ...
           'pFlow'], 'FontSize', 20)
    xlabel('time', 'FontSize', 20) % x-axis label
    ylabel('pFlow', 'FontSize', 20) % y-axis label
    h = legend('Numerical','Exact');
    set(h, 'FontSize', 16);

    str = sprintf('%s/LNFlow.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    btST = 0;           
    btET = 1.0;               
    rebTS = 0.01;


    simulator = ModelEvolver();
    simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                     model);

    F = size(simData, 1);
    m = size(simData, 2);

    maxRet = 0.02;
    maxDrawDown = -0.1;    
    maxLR = 2.0;
    constr = Constraint(maxRet, maxDrawDown, maxLR, false);
    
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
    
    plot(t, phiMat(1,:), 'Color', 'red', 'linewidth', 2);
    hold on;
    plot(t, exactStrategy(1,:), 'Color', 'blue', 'linewidth', 2);
    
    title('LN numerical and exact strategies', 'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20) % x-axis label
    ylabel('phi*', 'FontSize', 20) % y-axis label    
    h =legend('Numerical','Exact');
    set(h, 'FontSize', 16);
    
    axes('Position',[.17 .2 .2 .2])
    box on
    
    diff = phiMat(1,:) - exactStrategy;
    plot(t, diff, 'linewidth', 2);
    
    title('Difference', 'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20) % x-axis label
    ylabel('difference', 'FontSize', 20) % y-axis label    

    str = sprintf('%s/LNStrategy.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    % Ploting strategy and stock price
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

    str = sprintf('%s/LNPnPos.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    % Plotting PnL
    figure;
    plot(t, wVec(1,:), 'Color', 'red', 'linewidth', 2);    
    title('Cumulative PnL, LN', 'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20) % x-axis label
    ylabel('PnL', 'FontSize', 20) % y-axis label    

    str = sprintf('%s/LNPnL.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    toc
end

function runMR(outdir)

    tic
    
    display('Running MeanReverting');
    w0 = 100000;
    utilityType = 'HARA';
    turnedOnConsumption = false;
    numCores = 1;
    modelParam.modelType = 'MeanReverting';
    modelParam.mu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5];
    modelParam.vol = [0.1; 0.16; 0.3; 0.52; 0.14; 0.5; 1.0; 0.3; 0.5; 0.8];
    modelParam.lambda = [2.10; 1.32; 1.10; 1.24; 1.56; 0.6; 1.9; 2.3; 1.05; 0.8];

    F = length(modelParam.mu);

    % Generate correlation matrix
    tmp = randn(F, F);
    C = tmp * tmp';
    D = zeros(F, F);
    for i = 1:F
        C(i, i) = 2 * C(i, i);
        D(i,i) = sqrt(C(i,i));
    end
    corrMatr = inv(D) * C * inv(D)
    save('MeanRevertingcorr.mat', 'corrMatr');
    
    rankofCorr = rank(corrMatr)
    eigenV = eig(corrMatr)
    inverseCorr = inv(corrMatr)
    
    turnedOnConsumption = false;

    gamma = 10.0;
    xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1];
    pT = zeros(length(xCurr),1);
        
    tCurr = 0;
    T = 1.0;
    timeStep = 0.01;
    tol = 1e-6;
    

    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);    
    utiCalc = UtilityCalculator(gamma, utilityType);
    hamSys = HamiltonianSystem(portCalc, utiCalc);       
    wkbSolver = WKBHierarchySolver(hamSys, numCores);

    % Ploting Mean Reverting numerical flow
    xT = wkbSolver.solveForTermX(xCurr, tCurr, T, timeStep, tol);    
    [xFlow, pFlow] = wkbSolver.generateLfFlow(xT, pT, tCurr, T, timeStep, tol);
    
    % Ploting Mean Reverting exact flow
    t = tCurr:timeStep:T;    
    LAMBDA = diag(modelParam.lambda);
    cov = portCalc.instCov(xCurr);
    
    % changed due to the notation change, not the original ones
    A = [-hamSys.oneOverGamma * LAMBDA,  cov; ...
         -hamSys.kappa * hamSys.oneOverGamma * LAMBDA * inv(cov) * LAMBDA, hamSys.oneOverGamma ...
         * LAMBDA];
    
    m = [hamSys.oneOverGamma * LAMBDA * modelParam.mu; ...
         hamSys.kappa * hamSys.oneOverGamma * LAMBDA * inv(cov) * LAMBDA * modelParam.mu];
    
    [xFlowExact, pFlowExact] = calcExactMRFlow(t, A, T, xT, m);

    display('Plotting mean reverting exact xFlow and leapfrog flow on the same graph. They should be on top of each other') 
    
    F = length(xCurr);
    numAsset = 2;
    figure;        
    for i = 1:numAsset
        subplot(numAsset,2,2*i-1);
        plot(t, xFlow(i,:), 'Color', 'red', 'linewidth', 2);
        hold on;
        plot(t, xFlowExact(i,:), 'Color', 'blue', 'linewidth', 2);
        
        str1 = sprintf('Asset %d OU numerical xFlow and exact xFlow', i);
        title(str1, 'FontSize', 20);
        
        xlabel('time', 'FontSize', 20) % x-axis label
        ylabel('xFlow', 'FontSize', 20) % y-axis label    
        h = legend('Numerical','Exact');
        set(h, 'FontSize', 16);
        
        subplot(numAsset,2,2*i);
        plot(t, pFlow(i,:), 'Color', 'red', 'linewidth', 2);
        hold on;
        plot(t, pFlowExact(i,:), 'Color', 'blue', 'linewidth', 2);
        
        str2 = sprintf('Asset %d OU numerical pFlow and exact pFlow', i);
        title(str2, 'FontSize', 20);
        xlabel('time', 'FontSize', 20) % x-axis label
        ylabel('pFlow', 'FontSize', 20) % y-axis label        
        h = legend('Numerical', 'Exact');
        set(h, 'FontSize', 16);
    end
    
    str = sprintf('%s/MRFlow.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);

    btST = 0;           
    btET = 1.0;               
    rebTS = 0.02;

    simulator = ModelEvolver();
    simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                     model);

    maxRet = 0.2;
    maxDrawDown = -0.1;    
    maxLR = 2.0;
    constr = Constraint(maxRet, maxDrawDown, maxLR, true);

    bte = BtEngine(btST, btET, rebTS, constr);

    [wVec, phiMat, cVec] = bte.runBackTest(simData, wkbSolver, w0, ...
                                           turnedOnConsumption);

    t = btST:rebTS:btET;    
    
    % ploting strategy and simulated stock price
    figure;
    for i=1:numAsset
        subplot(numAsset,1, i)
        pos = [phiMat(i,:), 0];

        [ax, h1, h2] = plotyy(t, simData(i,:), t, pos, 'plot', 'stairs');    
        set(h1, 'linewidth', 2, 'color', 'red');
        set(h2,  'linewidth', 2, 'color', 'blue');
        
        str = sprintf(['Asset %d simulated stock price and wkb ' ...
                       'approximate position'],i);
        title(str,'FontSize', 20);
        
        xlabel('rebalance time', 'FontSize', 20) % x-axis label
        ylabel(ax(1), 'simulated stock price($)', 'color', 'red', 'FontSize', 20);
        ylabel(ax(2), 'position', 'color', 'blue', 'FontSize', 20);
        
        h = legend('simulated stock price','wkb strategy', 'Location', 'southeast');
        set(h, 'FontSize', 16);

    end

    str = sprintf('%s/MRPnPos.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    % Ploting PnL
    figure;
    plot(t, wVec(1,:), 'Color', 'red', 'linewidth', 2);    
    title('Cumulative PnL, OU', 'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20); 
    ylabel('PnL', 'FontSize', 20);  

    str = sprintf('%s/MRPnL.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    toc

end


function runCIR(outdir)

    tic
    
    display('Running CIR');
    w0 = 100000;
    utilityType = 'HARA';
    turnedOnConsumption = false;
    numCores = 1;
    modelParam.modelType = 'CIR';
    modelParam.mu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5];
    modelParam.vol = [0.1; 0.16; 0.3; 0.52; 0.14; 0.5; 1.0; 0.3; 0.5; 0.8];
    modelParam.vol = modelParam.vol ./ sqrt(modelParam.mu); 
    modelParam.lambda = [2.10; 1.32; 1.10; 1.24; 1.56; 0.6; 1.9; 2.3; 1.05; 0.8];
    xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1];
    
    F = length(modelParam.mu);

    tmp = randn(F, F);
    C = tmp * tmp';
    D = zeros(F, F);
    for i = 1:F
        C(i, i) = 2 * C(i, i);
        D(i,i) = sqrt(C(i,i));
    end
    
    corrMatr = inv(D) * C * inv(D)
    save('CIRcorr.mat', 'corrMatr');    
    rankofCorr = rank(corrMatr)
    eigenV = eig(corrMatr)
    inverseCorr = inv(corrMatr)
    
    turnedOnConsumption = false;

    gamma = 10.0;
    pT = zeros(length(xCurr),1);
        
    tCurr = 0;
    T = 1.0;
    timeStep = 0.01;
    tol = 1e-6;
    

    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);    
    utiCalc = UtilityCalculator(gamma, utilityType);
    hamSys = HamiltonianSystem(portCalc, utiCalc);       
    wkbSolver = WKBHierarchySolver(hamSys, numCores);

    % Ploting CIR numerical flow
    xT = wkbSolver.solveForTermX(xCurr, tCurr, T, timeStep, tol);    
    [xFlow, pFlow] = wkbSolver.generateLfFlow(xT, pT, tCurr, T, timeStep, tol);
    
    display('Plotting CIR leapfrog flow on the same graph.'); 
    
    F = length(xCurr);
    numAsset = 2;
    t = tCurr:timeStep:T; 
    figure;        
    for i = 1:numAsset
        subplot(numAsset,2,2*i-1);
        plot(t, xFlow(i,:), 'Color', 'red', 'linewidth', 2);
        
        str1 = sprintf('Asset %d CIR numerical xFlow', i);
        title(str1, 'FontSize', 20);
        
        xlabel('time', 'FontSize', 20) % x-axis label
        ylabel('xFlow', 'FontSize', 20) % y-axis label    
        
        subplot(numAsset,2,2*i);
        plot(t, pFlow(i,:), 'Color', 'red', 'linewidth', 2);
        
        str2 = sprintf('Asset %d CIR numericla pFlow', i);
        title(str2, 'FontSize', 20);
        xlabel('time', 'FontSize', 20) % x-axis label
        ylabel('pFlow', 'FontSize', 20) % y-axis label        
    end
    
    str = sprintf('%s/CIRFlow.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);

    btST = 0;           
    btET = 1.0;               
    rebTS = 0.02;

    simulator = ModelEvolver();
    simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                     model);

    maxRet = 0.2;
    maxDrawDown = -0.1;    
    maxLR = 2.0;
    constr = Constraint(maxRet, maxDrawDown, maxLR, true);
    
    bte = BtEngine(btST, btET, rebTS, constr);
    simData
    [wVec, phiMat, cVec] = bte.runBackTest(simData, wkbSolver, w0, ...
                                           turnedOnConsumption);

    t = btST:rebTS:btET;    
    
    % ploting strategy and simulated stock price
    figure;
    for i=1:numAsset
        subplot(numAsset,1, i)
        pos = [phiMat(i,:), 0];

        [ax, h1, h2] = plotyy(t, simData(i,:), t, pos, 'plot', 'stairs');    
        set(h1, 'linewidth', 2, 'color', 'red');
        set(h2,  'linewidth', 2, 'color', 'blue');
        
        str = sprintf(['Asset %d simulated stock price and wkb ' ...
                       'approximate position'],i);
        title(str,'FontSize', 20);
        
        xlabel('rebalance time', 'FontSize', 20) % x-axis label
        ylabel(ax(1), 'simulated stock price($)', 'color', 'red', 'FontSize', 20);
        ylabel(ax(2), 'position', 'color', 'blue', 'FontSize', 20);
        
        h = legend('simulated stock price','wkb strategy', 'Location', 'southeast');
        set(h, 'FontSize', 16);

    end

    str = sprintf('%s/CIRPnPos.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);
    saveas(gcf, str);
    
    % Ploting PnL
    figure;
    plot(t, wVec(1,:), 'Color', 'red', 'linewidth', 2);    
    title('Cumulative PnL, CIR', 'FontSize', 20);
    xlabel('rebalance time', 'FontSize', 20); 
    ylabel('PnL', 'FontSize', 20);  

    str = sprintf('%s/CIRPnL.png', outdir);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                       20]);    
    saveas(gcf, str);
    
    toc

end


function [xFlow, pFlow] = calcExactMRFlow(sVec, A, T, y, m)

    numPoint = length(sVec);
    
    
    n = length(y);
    zerosVec = zeros(n, 1);
    termVal = [y; zerosVec];
    
    res = zeros(2*n, numPoint);
    
    for i = 1:numPoint
        
        s = sVec(i);
        term1 = expm(-(T - s) * A);
        term2 = (termVal + inv(A) * m);
        term3 = inv(A) * m;
        res(:,i) = term1 * term2 - term3;
        
    end
    
    xFlow = res(1:n, :);
    pFlow = res((n+1):2*n, :);
    
end


function runMonteCarlo(numPaths, modelParam, xCurr)

    str = sprintf('Monte Carlo back test, %s', modelParam.modelType);
    display(str);
    
    w0 = 100000;
    utilityType = 'HARA';
    turnedOnConsumption = false;
    numCores = 1;
    gamma = 10;
    
    btST = 0;           
    btET = 1.0;               
    rebTS = 0.02;

    maxRet = 0.2;
    maxDrawDown = -0.1;    
    maxLR = 2.0;

    % number of asset to plot the simulated stock price and pos
    numAsset = 2;
    
    corrMatrFile = strcat(modelParam.modelType, 'corr.mat');
    load(corrMatrFile);
    
    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);    
    utiCalc = UtilityCalculator(gamma, utilityType);
    hamSys = HamiltonianSystem(portCalc, utiCalc);       
    wkbSolver = WKBHierarchySolver(hamSys, numCores);

    wTVec = zeros(1, numPaths);
    UTVec = zeros(1, numPaths);
    
    for path = 1:numPaths
        
        tic
        str = sprintf('Running model %s path %d', modelParam.modelType, ...
                      path);
        display(str);
        
        simulator = ModelEvolver();
        simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                        model);

        constr = Constraint(maxRet, maxDrawDown, maxLR, true);
        
        bte = BtEngine(btST, btET, rebTS, constr);
        
        [wVec, phiMat, cVec] = bte.runBackTest(simData, wkbSolver, w0, ...
                                               turnedOnConsumption);
        
        wTVec(path) = wVec(end);
        UTVec(path) = utiCalc.U( wVec(end) );

        % Start to plot the back-test result
        t = btST:rebTS:btET;    
        
        % ploting strategy and simulated stock price
        figure;
        for i=1:numAsset
            subplot(numAsset,1, i)
            pos = [phiMat(i,:), 0];

            [ax, h1, h2] = plotyy(t, simData(i,:), t, pos, 'plot', 'stairs');    
            set(h1, 'linewidth', 2, 'color', 'red');
            set(h2,  'linewidth', 2, 'color', 'blue');
            
            str = sprintf(['Asset %d simulated stock price and wkb ' ...
                           'approximate position'],i);
            title(str,'FontSize', 20);
            
            xlabel('rebalance time', 'FontSize', 20) % x-axis label
            ylabel(ax(1), 'simulated stock price($)', 'color', 'red', 'FontSize', 20);
            ylabel(ax(2), 'position', 'color', 'blue', 'FontSize', 20);
            
            h = legend('simulated stock price','wkb strategy', 'Location', 'southeast');
            set(h, 'FontSize', 16);

        end

        str = sprintf('MC/%sPnPosPath%d.png', modelParam.modelType, ...
                      path);
        
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                            20]);
        saveas(gcf, str);
        
        % Ploting PnL
        figure;
        plot(t, wVec(1,:), 'Color', 'red', 'linewidth', 2);
        
        str = sprintf('Cumulative PnL %s', modelParam.modelType);
        title(str, 'FontSize', 20);
        xlabel('rebalance time', 'FontSize', 20); 
        ylabel('PnL', 'FontSize', 20);  

        str = sprintf('MC/%sPnLPath%d.png', modelParam.modelType, path);
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 ...
                            20]);    
        saveas(gcf, str);
        
        toc
    
    end
    
    averageWT = sum(wTVec) / numPaths
    averageUT = sum(UTVec) / numPaths
    
    hist(wTVec)
    str = sprintf('MC WT distribution %s', modelParam.modelType);
    title(str, 'FontSize', 20);
    xlabel('WT($)', 'FontSize', 20); 
    ylabel('Frequency', 'FontSize', 20);  

    str = sprintf('MC/%sWTHist.png', modelParam.modelType);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 20]);    
    saveas(gcf, str);
    
    hist(UTVec)
    str = sprintf('MC UT distribution %s', modelParam.modelType);
    title(str, 'FontSize', 20);
    xlabel('UT($)', 'FontSize', 20); 
    ylabel('Frequency', 'FontSize', 20);  

    str = sprintf('MC/%sUTHist.png', modelParam.modelType);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 20]);    
    saveas(gcf, str);
    
    str = sprintf('MC/%sWTVec.mat', modelParam.modelType);
    save(str, 'wTVec');

    str = sprintf('MC/%sUTVec.mat', modelParam.modelType);    
    save(str, 'UTVec');

end
