clear;clc;

% Global variables

global isMeanReverting;
isMeanReverting = false   % if it is true, it will run a mean reverting
                         % model, otherwise, lognomal model 

% Initialize model parameters
if isMeanReverting

  % MEAN REVERTING model parameters. 
  global MRMu;
  MRMu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5]      % Dim: same as x, COLUMN
  global lambda;
  lambda = [21.0; 13.2; 11.0; 12.4; 15.6; 6; 19; 23; 10.5; 8]  % Dim: same as x, COLUMN
  global nVols;
  nVols = [0.01; 0.016; 0.03; 0.052; 0.014; 0.05; 0.10; 0.03; 0.05; 0.08]  % We are assuming each asset have only one
                                                                 % Brownian motion, so this will be used to
                                                                 % generated a diagonal matrix b. For simplicity,
                                                                 % we are using n=p at this point. 

else

  % LOGNORMAL model parameters
  global lnMu;
  lnMu = [0.8; 1.3; 0.9]    % Dim: same as x, COLUMN 
  global lnVols;
  lnVols = [0.33; 0.16; 0.45] % We are assuming each asset have only one 
                        % Brownian motion, so this will be used to
                        % generated a diagonal matrix b. For simplicity,
                        % we are using n=p at this point.  
end

global corrMatr;     
corrMatr = eye(3)  % correlation matrix between dZs, assuming they are
                   % independent. Dim: p*p
gamma = 10.0;
global oneOverGamma;
oneOverGamma = 1 /gamma
global kappa;
kappa = oneOverGamma - 1

xCurr = [1.0; 0.8; 1.5]
%xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1]  % asset price vector
tCurr = 0           
T = 10               
timeStep = 0.01
tol = 0.0001
maxIter = 1000000000   % make it very large so the iteration conveges

% Test [phi] = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter) function
phi = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter)

phi_exact = invInstCov(xCurr) * driftV(xCurr)

diff = norm(phi - phi_exact) 

xT = solveForTermX(xCurr, tCurr, T, timeStep, tol, maxIter)
% expectedProfit = phi' * (xT - xCurr)


[xFlow, pFlow] = generateLfFlow(xT, [0;0;0], tCurr, T, timeStep, tol, maxIter);
t = tCurr:timeStep:T;


x_exact_1 = xT(1) * exp(-lnMu(1) / gamma * (T - t));
x_exact_2 = xT(2) * exp(-lnMu(2) / gamma * (T - t));
x_exact_3 = xT(3) * exp(-lnMu(3) / gamma * (T - t));

figure;
subplot(2,2,1);
plot(t, xFlow(1,:)); 
hold on;
plot(t, x_exact_1);

subplot(2,2,2);
plot(t, xFlow(2,:));
hold on;
plot(t, x_exact_2);

subplot(2,2,3);
plot(t, xFlow(3,:));
hold on;
plot(t, x_exact_3);

pause;
