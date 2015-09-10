% This is a main script for funtions test.
% Author: Zhenfeng Liang, Wei Liu
% Contact: Zhenfeng.Jose.Liang@gmail.com
% All copyrights(C) reserved.

% The following test is assuming two dimensions asset price.
x = [1, 1];  % asset price vector
corrMatr = [1,0.5;0.5,1];

% Test driftV function.
drofx = driftV(x)

% Test diffV function
dofx = diffV(x)

% Test instCov function
#covx = instCov(x, corrMatr) 

% Test diffDer function
dDofx = diffDer(x, 1)

% Test driftDer function
drDofx = driftDer(x,1)

% Test driftDer2 function
drD2ofx = driftDer2(x, 1, 2)

% Test diffDer2 function
dD2ofx = diffDer2(x, 1, 2)
