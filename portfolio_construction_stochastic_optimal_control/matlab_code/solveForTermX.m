function [xT] = solveForTermX(x, t, T, timeStep, tol, maxIter)
  % Invert the Hamiltonian flow (zero terminal momentum)
  % Input: x, initial asset price vector, t, starting time, T, terminal time, timeStep: time step for flow, tol: tolerance, maxIter: maximum iteration
  % Output: xT, terminal value y
  % Note: Original R code messed up dimension of the Flow when it tries to use Newton method to solve it.

  F = length(x);

  % setup the system of approximate variational equations
  zeroVec = zeros(F, 1);
  termVal = [eye(F); zeros(F)];

  
  xT = x;    
  err = 1.0;
  ctr = 1;    
  while err >= tol && ctr <= maxIter
        
    [xFlow, pFlow] = generateLfFlow(xT, zeroVec, t, T, timeStep, tol, maxIter);

    Q = hamDxp(xT, zeroVec);
    R = hamDpp(xT, zeroVec);
    U = hamDxx(xT, zeroVec);
    
    apprVarSol = expm((t-T) * [Q, R; -U, -Q]) * termVal;
    invDF = inv(apprVarSol(1:F, :));

    z = xT - invDF * (xFlow(1:F, 1) - x);         % Note: original R code messed up dimensions here!!!
    err = norm(z - xT);
    xT = z;
    ctr = ctr + 1;

  end

  display('DONE INVERTING THE FLOW')
%  display(xT)
  
end
