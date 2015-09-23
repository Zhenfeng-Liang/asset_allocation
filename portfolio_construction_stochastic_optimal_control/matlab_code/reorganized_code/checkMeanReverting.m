function checkMeanReverting()

    modelParam.modelType = 'MeanReverting';
    modelParam.mu = [0.4; 1.3; 2.2; 3.5];
    modelParam.vol = [0.01; 0.016; 0.03; 0.052];
    modelParam.lambda = [21.0; 13.2; 11.0; 12.4];
    corrMatr = eye(4);

    gamma = 10.0;
    xCurr = [0.8; 0.8; 2; 4];
    tCurr = 0;           
    T = 1;               
    timeStep = 0.01;
    tol = 0.0001;
    
    model = Model(modelParam);
    portCalc = PortfolioCalculator(model, corrMatr);
    hamSys = HamiltonianSystem(portCalc, gamma);
    wkbSolver = WKBHierarchySolver(hamSys);
    
    y = wkbSolver.solveForTermX(xCurr, tCurr, T, timeStep, tol);

    t = tCurr:timeStep:T;
    
    LAMBDA = diag(modelParam.lambda);
    cov = portCalc.instCov(xCurr);
    
    
    A = [-hamSys.oneOverGamma * LAMBDA, hamSys.oneOverGamma * cov; ...
         -hamSys.kappa * LAMBDA * inv(cov) * LAMBDA, hamSys.oneOverGamma ...
         * LAMBDA];
    
    m = [hamSys.oneOverGamma * LAMBDA * modelParam.mu; ...
         hamSys.kappa * LAMBDA * inv(cov) * LAMBDA * modelParam.mu];
    
    [xFlowExact, pFlowExact] = calcExactMRFlow(t, A, T, y, m);

    [xFlow, pFlow] = wkbSolver.generateLfFlow(y, zeros(4,1), ...
                                              tCurr, T, timeStep, tol);
    
    
    display('Plotting mean reverting exact xFlow and leapfrog flow on the same graph. They should be on top of each other') 
    
    figure;
    
    subplot(2,2,1);
    plot(t, xFlow(1,:), 'Color', 'red');
    hold on;
    plot(t, xFlowExact(1,:), 'Color', 'blue');
    
    title('Asset 1: mean reverting wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label    
    
    subplot(2,2,2);
    plot(t, xFlow(2,:), 'Color', 'red');
    hold on;
    plot(t, xFlowExact(2,:), 'Color', 'blue');
    
    title('Asset 2: mean reverting wkb approximate xFlow and exact xFlow');
    xlabel('time') % x-axis label
    ylabel('xFlow') % y-axis label
    
    subplot(2,2,3);
    plot(t, pFlow(1,:), 'Color', 'red');
    hold on;
    plot(t, pFlowExact(1,:), 'Color', 'blue');

    title('Asset 1: mean reverting wkb approximate pFlow and exact pFlow');
    xlabel('time') % x-axis label
    ylabel('pFlow') % y-axis label

    
    subplot(2,2,4);
    plot(t, pFlow(2,:), 'Color', 'red');
    hold on;
    plot(t, pFlowExact(2,:), 'Color', 'blue');

    title('Asset 2: mean reverting wkb approximate pFlow and exact pFlow');
    xlabel('time') % x-axis label
    ylabel('pFlow') % y-axis label

    
    display(['Checking strategies with and without S1, when correlation ' ...
             'matrix is an identity matrix they should be the same']);
    
    independentWithS1 = 1.0 / gamma * wkbSolver.optimalControlStrategy(xCurr, ...
                                                      tCurr, T, timeStep, tol)
    wkbSolver.includeS1 = false;
    
    independentWithoutS1 = 1.0 / gamma * ...
        wkbSolver.optimalControlStrategy(xCurr, tCurr, T, timeStep, tol)
    
    diff = independentWithoutS1 - independentWithS1
    diffNorm = norm(diff)

    display('Finished checking strategy for independent case.')


    
    display(['Checking strategies with and without S1, when correlation ' ...
             'matrix is not an identity matrix. We can see the ' ...
             'dependent assets are no longer the same but the ' ...
             'independent ones are which is inconsistent with the ' ...
             'paper argument']);
    
    corrMatr(1,2) = 0.5;
    corrMatr(2,1) = 0.5;
    corrMatr(1,3) = -0.1;
    corrMatr(3,1) = -0.1;
    corrMatr(2,3) = 0.3;
    corrMatr(3,2) = 0.3;
    
    portCalc = PortfolioCalculator(model, corrMatr);
    hamSys = HamiltonianSystem(portCalc, gamma);
    wkbSolver = WKBHierarchySolver(hamSys);
    
    dependentWithS1 = 1.0 / gamma * wkbSolver.optimalControlStrategy(xCurr, ...
                                                      tCurr, T, timeStep, tol)
    wkbSolver.includeS1 = false;
    
    dependentWithoutS1 = 1.0 / gamma * ...
        wkbSolver.optimalControlStrategy(xCurr, tCurr, T, timeStep, tol)
    
    diff = dependentWithoutS1 - dependentWithS1
    diffNorm = norm(diff)

    display('Finished checking strategy for dependent case.')

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
