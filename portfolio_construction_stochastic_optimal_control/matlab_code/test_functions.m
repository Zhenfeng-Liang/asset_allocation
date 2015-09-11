% This is a main script for funtions test.
% Author: Zhenfeng Liang, Wei Liu
% Contact: Zhenfeng.Jose.Liang@gmail.com
% All copyrights(C) reserved.

% The following test is assuming two dimensions asset price.

# Global variables
global corrMatr;
corrMatr = [1,0.5;0.5,1];
gamma = 0.5;
global oneOverGamma = 1 /gamma;
global kappa = oneOverGamma - 1;

x = [1; 1]  # asset price vector
p = [1;1]   # Have to be n x 1 dim

% Test driftV function.
drofx = driftV(x)

% Test diffV function
dofx = diffV(x)

% Test instCov function
covx = instCov(x) 

% Test diffDer function
dDofx = diffDer(x, 1)

% Test driftDer function
drDofx = driftDer(x,1)

% Test driftDer2 function
drD2ofx = driftDer2(x, 1, 2)

% Test diffDer2 function
dD2ofx = diffDer2(x, 1, 2)

% Test invInstCov function
iCovx = invInstCov(x)

% Test instCovDer function
covxDi = instCovDer(x, 1)

% Test invInstCovDer function
iCovxDi = invInstCovDer(x, 1)

% Test instCovDer2 function
covxDij = instCovDer2(x, 1, 2)

% Test invInstCovDer2 function
iCovxDij = invInstCovDer2(x, 1, 2)

% Test lagr(x, p) function
lagrval = lagr(x, p)

% Test ham(x, p) function
hVal = ham(x, p)

% Test hamDx(x, p) function
hDx = hamDx(x, p)

% Test hamDp(x, p) function
hDp = hamDp(x, p)

% Test hamDxx(x, p) function
hDxx = hamDxx(x, p)

% Test hamDxp(x, p) function
hDxp = hamDxp(x, p)

% Test hamDpp(x, p) function
hDpp = hamDpp(x, p)
