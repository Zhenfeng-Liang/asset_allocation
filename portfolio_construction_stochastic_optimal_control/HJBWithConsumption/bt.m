close all;
clc;
clear;

tic
w0 = 1e5;
utilityType = 'CRRA';
turnedOnConsumption = false;
numCores = 2;


%modelParam.modelType = 'CIR';
%modelParam.mu = [0.4; 0.6];
%modelParam.vol = [0.18; 0.16];
%modelParam.lambda = [2.0; 1.5];
%

%modelParam.modelType = 'LogNormal'
%modelParam.mu = [0.8; 1.3];
%modelParam.vol = [0.33; 0.16];

modelParam.modelType = 'MeanReverting';
modelParam.mu = [0.4; 0.6];
modelParam.vol = [0.15; 0.16];
modelParam.lambda = [2.0; 1.5];

corrMatr = eye(2.0);

gamma = 10.0;
tol = 1e-6;

model = Model(modelParam);
portCalc = PortfolioCalculator(model, corrMatr);

utiCalc = UtilityCalculator(gamma, utilityType);

hamSys = HamiltonianSystem(portCalc, utiCalc);

wkbSolver = WKBHierarchySolver(hamSys, numCores);


% back test parameters
btST = 0;           
btET = 0.5;               
rebTS = 0.01;
turnedOnConsumtion = false;

xCurr = [0.8; 0.8];

simulator = ModelEvolver();
simData = simulator.EvolveEuler(xCurr, btST, btET, rebTS, corrMatr, ...
                                     model)
            
%histData = [0.8, 0.6, 0.7, 0.5, 0.4, 0.3, 0.5;
%            0.8, 0.9, 0.85, 0.7, 1.0, 1.1, 1.3];


constr = Constraint(0.3, -0.2, false);
bte = BtEngine(btST, btET, rebTS, constr);

[wVec, phiMat, cVec] = bte.runBackTest(simData, wkbSolver, w0, ...
                                       turnedOnConsumtion)


toc
