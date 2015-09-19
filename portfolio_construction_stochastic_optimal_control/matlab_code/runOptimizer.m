clear;clc;

% Global variables

global isMeanReverting;
isMeanReverting = true   % if it is true, it will run a mean reverting
                         % model, otherwise, lognomal model 

% Initialize model parameters
if isMeanReverting

  % MEAN REVERTING model parameters. 
  global MRMu;
  MRMu = [0.4; 1.3; 2.2; 3.5; 1.2; 4.0; 5.5; 2.0; 1.0; 4.5]      % Dim: same as x, COLUMN
  global lambda;
  lambda = [21.0; 13.2; 11.0; 12.4; 15.6; 6; 19; 23; 10.5; 8]  % Dim: same as x, COLUMN
  global nVols;
  nVols = [0.1; 0.16; 0.3; 0.52; 0.14; 0.5; 1.0; 0.3; 0.5; 0.8]  % We are assuming each asset have only one
                                                                 % Brownian motion, so this will be used to
                                                                 % generated a diagonal matrix b. For simplicity,
                                                                 % we are using n=p at this point. 

else

  % LOGNORMAL model parameters
  global lnMu;
  lnMu = [0.4; 0.12]    % Dim: same as x, COLUMN 
  global lnVols;
  lnVols = [0.33; 0.26] % We are assuming each asset have only one 
                        % Brownian motion, so this will be used to
                        % generated a diagonal matrix b. For simplicity,
                        % we are using n=p at this point.  
end

global corrMatr;     
corrMatr = eye(10)  % correlation matrix between dZs, assuming they are
                   % independent. Dim: p*p
gamma = 10;
global oneOverGamma;
oneOverGamma = 1 /gamma
global kappa;
kappa = oneOverGamma - 1


xCurr = [0.8; 0.8; 2; 4; 1; 3; 6; 1.5; 2.4; 5.1]  % asset price vector
tCurr = 0           
T = 1               
timeStep = 0.01
tol = 0.0001
maxIter = 1000000000   % make it very large so the iteration conveges

% Test [phi] = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter) function
phi = optimalControlStrategy(xCurr, tCurr, T, timeStep, tol, maxIter)



