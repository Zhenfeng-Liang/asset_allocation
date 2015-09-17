clear;clc;
% Global variables
global corrMatr;
corrMatr = [1.0,0; 0,1.0]
gamma = 10;
global oneOverGamma;
oneOverGamma = 1 /gamma
global kappa;
kappa = oneOverGamma - 1

x = [0.8; 0.8]  % asset price vector

t = 0
T = 1
timeStep = 0.01
tol = 0.0001
maxIter = 1000000000

% Test [phi] = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter) function
phi = optimalControlStrategy(x, t, T, timeStep, tol, maxIter)



