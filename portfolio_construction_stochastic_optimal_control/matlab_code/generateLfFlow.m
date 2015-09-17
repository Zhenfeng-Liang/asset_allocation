function [xFlow, pFlow] = generateLfFlow(xT, pT, t, T, timeStep, tol, maxIter)
  % Input: xT, pT, terminal value vector of asset price and momentum vector, from t to T, supposed timeStep, tol: error tolerance for leapfrog implicit equation, maxIter: maximum iteration for implicit iteration
  % Output: x(t) and p(t) flow map given xT and pT

  numSteps = ceil((T - t) / timeStep);
  h_over_2 = 0.5 * (T - t) / numSteps;  % Note: original R code doesn't have t here. I suspect that the code assume t = 0
  
  numFac = length(xT);

  xFlow = zeros(numFac, numSteps + 1); % Note: original R code use dimension numSteps here, in that case, what it returns is the value at (t + h)
  pFlow = zeros(numFac, numSteps + 1); % Note: original R code use dimension numSteps here, in that case, what it returns is the value at (t + h)
  
  x = xT;
  p = pT;
  for i = (numSteps + 1):-1:2  % Note: original R code starts from numSteps

    xFlow(:,i) = x;
    pFlow(:,i) = p;
    [x, p] = oneTimeStep(x, p, h_over_2, tol, maxIter);

  end

  xFlow(:,1) = x;
  pFlow(:,1) = p;

end



function [x_minus_1, p_minus_1] = oneTimeStep(x, p, h_over_2, tol, maxIter)
  % Input: x, p: vector value at this point, h_over_2: half increament of time step, tol: error tolerance for leapfrog implicit equation, maxIter: maximum iteration for implicit iteration
  % Output: one time backward x and p vector in the flow map

  p_minus_half_old = p;       % Initial guess
  p_minus_half_new = p;       % Declare the vector
  x_minus_1_old = x;          % Initial guess
  x_minus_1_new = x;          % Declare the vector
  
  err = 1.0;
  iter = 0;
  while (err > tol) && (iter <= maxIter)
        
    % The following lines uses Newton iteration
    f = p + h_over_2 * hamDx(x, p_minus_half_old) - p_minus_half_old;
    df = h_over_2 * hamDxp(x, p_minus_half_old) - eye(length(x));
    p_minus_half_new = p_minus_half_old - inv(df) * f;

    %p_minus_half_new = p + h_over_2 * hamDx(x, p_minus_half_old);          % this line is the original R using normal iteration

    err = norm(p_minus_half_new - p_minus_half_old);
    p_minus_half_old = p_minus_half_new;
    iter = iter + 1;
  end
  

  err = 1.0;
  iter = 0;
  while (err > tol) && (iter <= maxIter)
        
    % The following lines uses Newton iteration
    f = x_minus_1_old + h_over_2 * (hamDp(x, p_minus_half_new) + hamDp(x_minus_1_old, p_minus_half_new)) - x;
    df = h_over_2 * hamDxp(x_minus_1_old, p_minus_half_new) + eye(length(x));
    x_minus_1_new = x_minus_1_old - inv(df) * f;

    %x_minus_1_new = x - h_over_2 * (hamDp(x, p_minus_half_new) + hamDp(x_minus_1_old, p_minus_half_new));   % this line is the original R using normal iteration

    err = norm(x_minus_1_new - x_minus_1_old);
    x_minus_1_old = x_minus_1_new;
    iter = iter + 1;
  end
  
  x_minus_1 = x_minus_1_new;
  p_minus_1 = p_minus_half_new + h_over_2 * hamDx(x_minus_1_new, p_minus_half_new);
  
end



