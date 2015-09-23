close all;
clc;
clear;


% This is checking lognormal model with its exact solution.
checkLogNormal();

% This is checking mean reverting model with its exact solution
checkMeanReverting();


% Optimize the mean reverting model
modelParam.modelType = 'MeanReverting';
modelParam.mu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5];
modelParam.vol = [0.01; 0.016; 0.03; 0.052; 0.014; 0.05; 0.10; 0.03; 0.05; 0.08];
modelParam.lambda = [21.0; 13.2; 11.0; 12.4; 15.6; 6; 19; 23; 10.5; 8];
corrMatr = eye(10);
gamma = 10.0;
xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1];
tCurr = 0;           
T = 1;               
timeStep = 0.01;
tol = 0.0001;

strategy = wkbOptimizer(modelParam, corrMatr, gamma, xCurr, tCurr, T, timeStep, tol)


