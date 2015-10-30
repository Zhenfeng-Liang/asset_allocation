close all;
clc;
clear;

tic
w0 = 1e5;
utilityType = 'CRRA';
turnedOnConsumption = true;
numCores = 2;
%modelParam.modelType = 'MeanReverting';
%modelParam.modelType = 'LogNormal'

modelParam.modelType = 'CIR';
modelParam.mu = [0.4; 0.6];
modelParam.vol = [0.18; 0.16];
modelParam.lambda = [2.0; 1.5];

%modelParam.mu = [0.8; 1.3];
%modelParam.vol = [0.33; 0.16];
%modelParam.mu = [0.4; 1.3];
%modelParam.vol = [0.01; 0.016];
%modelParam.lambda = [0.210; 0.132];

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
rebTS = 0.1;
turnedOnConsumtion = true;

histData = [0.8, 0.6, 0.7, 0.5, 0.4, 0.3, 0.5;
            0.8, 0.9, 0.85, 0.7, 1.0, 1.1, 1.3];


bte = btEngine(btST, btET, rebTS);
[wVec, phiMat, cVec] = bte.runBackTest(histData, wkbSolver, w0, ...
                                       turnedOnConsumtion)

toc
