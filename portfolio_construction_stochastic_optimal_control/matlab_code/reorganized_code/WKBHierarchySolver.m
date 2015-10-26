classdef WKBHierarchySolver
    % This class implements the main part of WKB solver 
    
    
    properties
        hamSys;
        includeS1;
    end
    
    methods
    
        function obj = WKBHierarchySolver(hamSys)
            % hamSys: a hanSysCalc object
            
            obj.hamSys = hamSys;
            obj.includeS1 = true;
        end
        
        function [phi] = optimalControlStrategy(obj, xCurr, tCurr, T, timeStep, tol)
            % Input: xCurr: current asset price vector, tCurr: current
            % time, T: terminal time, timeStep: time interval between
            % points, tol: tolerance, maxIter: maximum iteration
            % number            
            % Output: phi: asset allocation weight vector
            % Note: make sure (T - tCurr) / timeStep is even, i.e. numStep
            % is odd. 
            
            phi = obj.hamSys.portCalc.invInstCov(xCurr) * obj.hamSys.portCalc.model.driftV(xCurr) ...
                + obj.calcNablaS(xCurr, tCurr, T, timeStep, tol);
            
        end

        function [xT] = solveForTermX(obj, x, t, T, timeStep, tol)
            % Invert the Hamiltonian flow (zero terminal momentum)
            % Input: x, initial asset price vector, t, starting time, T,
            % terminal time, timeStep: time step for flow, tol: tolerance,
            % maxIter: maximum iteration  
            % Output: xT, terminal value y
            % Note: Original R code messed up dimension of the Flow when it
            % tries to use Newton method to solve it. 
            
            F = length(x);
            
            % setup the system of approximate variational equations
            zeroVec = zeros(F, 1);
            termVal = [eye(F); zeros(F)];
            
            xT = x;
            err = 1.0;
            ctr = 1;
            while err >= tol
                
                [xFlow, pFlow] = obj.generateLfFlow(xT, zeroVec, t, T, timeStep, tol);
                
                Q = obj.hamSys.hamDxp(xT, zeroVec);
                R = obj.hamSys.hamDpp(xT, zeroVec);
                U = obj.hamSys.hamDxx(xT, zeroVec);
                
                apprVarSol = expm((t-T) * [Q, R; -U, -Q]) * termVal;
                invDF = inv(apprVarSol(1:F, :));
                
                z = xT - invDF * (xFlow(1:F, 1) - x);         % Note: original R code messed up dimensions here!!!
                err = norm(z - xT);
                xT = z;
                ctr = ctr + 1;
                
            end
            
            display('DONE INVERTING THE FLOW')
        end        
    
        
        function [S] = calcS(obj, x, t, T, timeStep, tol)
            % Input: x, initial asset price vector, t: starting time, T:
            % terminal time
            % Output: S: action from WKB approximation
            % Note: make sure (T - t) / timeStep is odd, i.e. numStep is
            % odd. Right now I didn't include S1 term.
            
            F = length(x);
            termVal = obj.solveForTermX(x, t, T, timeStep, tol);
            zeroVec = zeros(F, 1);
            
            [xPath, pPath] = obj.generateLfFlow(termVal, zeroVec, t, T, timeStep, tol);
            
            numSteps = size(xPath, 2) - 1;   % number of step equals to number of point minus one.
            
            dpDx = zeros(F, F, numSteps);
            for i = 1:(numSteps)
                
                one_over_dx = 1 ./ (xPath(:,i+1) - xPath(:,i));
                dp = pPath(:,i+1) - pPath(:,i);
                dpDx(:, :, i) = dp * one_over_dx';  % outer product
                
            end
            
            
            S0 = 0.0;
            S1 = 0.0;
            for i = 1:(size(xPath, 2)/2 - 1)   % Note the original R code is different here
                
                S0 = S0 + obj.lagr(xPath(:,2*i-1), pPath(:, 2*i-1)) ...
                    + 4 * obj.lagr(xPath(:, 2*i), pPath(:, 2*i)) ...
                    + obj.lagr(xPath(:, 2*i+1), pPath(:, 2*i+1));
                
                S1 = S1 + trace(obj.hamSys.portCalc.instCov(xPath(:,2*i-1)) * dpDx(:,:,2*i-1) ...
                    + 4 * obj.hamSys.portCalc.instCov(xPath(:,2*i)) * dpDx(:,:,2*i) ...
                    + obj.hamSys.portCalc.instCov(xPath(:,2*i+1)) * dpDx(:,:,2*i+1));
                
            end
            
            S0 = S0 * timeStep / 3
            S1 = S1 * timeStep / 6;   % S1 have another 1/2 coeff in front of
            % the integral
                        
            if obj.includeS1
                S = S0 + 1.0 * S1;
            else
                S = S0;
            end
            
        end
        
        
        function [nablaS] = calcNablaS(obj, x, t, T, timeStep, tol)
            % Input: x: initial asset price vector, t: starting time, T:
            % terminal time, timeStep: time interval between points, tol:
            % tolerance, maxIter: maximum iteration times.  
            % Output: gradient of S with respect to x vector
            % Note: make sure (T - t) / timeStep is odd.
            
            
            F = length(x);
            S = obj.calcS(x, t, T, timeStep, tol);
            nablaS = zeros(F, 1);
            eps = 1e-5;
            
            for i = 1:F
                xEps = x;
                xEps(i) = xEps(i) + eps;
                nablaS(i) = (obj.calcS(xEps, t, T, timeStep, tol) - S) / eps;
            end
            
        end

        function [val] = lagr(obj, x, p)
            % Lagrange function
            % Input: x, asset price vector, p vector,
            % Output: lagrange value
            % Note: In the R code, it uses term1 - term2, I am not sure if it is what it should be. I am suspecting that the original code ignore the minus sign before According to the paper, this should be term2 - term1, which is what I am using here. this need to be double checked later.
            
            a = obj.hamSys.portCalc.model.driftV(x);
            
            term1 = 0.5 * obj.hamSys.oneOverGamma * p' * obj.hamSys.portCalc.instCov(x) * p;
            term2 = 0.5 * obj.hamSys.kappa * a' * obj.hamSys.portCalc.invInstCov(x) * a;
            val = term2 - term1;
            
        end

        function [xFlow, pFlow] = generateLfFlow(obj, xT, pT, t, T, timeStep, tol)
            % Input: xT, pT, terminal value vector of asset price and
            % momentum vector, from t to T, supposed timeStep, tol: error
            % tolerance for leapfrog implicit equation, maxIter: maximum
            % iteration for implicit iteration   
            % Output: x(t) and p(t) flow map given xT and pT
            
            numSteps = ceil((T - t) / timeStep);
            
            % Note: original R code doesn't have t here. I suspect that the
            % code assume t = 0 
            h_over_2 = 0.5 * (T - t) / numSteps;  
            
            numFac = length(xT);
            
            xFlow = zeros(numFac, numSteps + 1); 
            pFlow = zeros(numFac, numSteps + 1); 
            
            x = xT;
            p = pT;
            for i = (numSteps + 1):-1:2  
                
                xFlow(:,i) = x;
                pFlow(:,i) = p;
                [x, p] = obj.oneTimeStep(x, p, h_over_2, tol);
                
            end
            
            xFlow(:,1) = x;
            pFlow(:,1) = p;
            
        end
        
        
        
        function [x_minus_1, p_minus_1] = oneTimeStep(obj, x, p, h_over_2, tol)
            % Input: x, p: vector value at this point, h_over_2: half
            % increament of time step, tol: error tolerance for leapfrog
            % implicit equation, maxIter: maximum iteration for implicit
            % iteration   
            % Output: one time backward x and p vector in the flow map
            
            p_minus_half_old = p;       % Initial guess
            p_minus_half_new = p;       % Declare the vector
            x_minus_1_old = x;          % Initial guess
            x_minus_1_new = x;          % Declare the vector
            
            err = 1.0;
            iter = 0; 
            while (err > tol)
                
                % The following lines uses Newton iteration
                f = p + h_over_2 * obj.hamSys.hamDx(x, p_minus_half_old) - p_minus_half_old;
                df = h_over_2 * obj.hamSys.hamDxp(x, p_minus_half_old) - eye(length(x));
                p_minus_half_new = p_minus_half_old - inv(df) * f;
                                
                err = norm(p_minus_half_new - p_minus_half_old);
                p_minus_half_old = p_minus_half_new;
                iter = iter + 1;
            end
            
            
            err = 1.0;
            iter = 0; 
            while (err > tol)
                
                % The following lines uses Newton iteration
                f = x_minus_1_old + h_over_2 * (obj.hamSys.hamDp(x, p_minus_half_new) + obj.hamSys.hamDp(x_minus_1_old, p_minus_half_new)) - x;
                df = h_over_2 * obj.hamSys.hamDxp(x_minus_1_old, p_minus_half_new) + eye(length(x));
                x_minus_1_new = x_minus_1_old - inv(df) * f;
                
                
                err = norm(x_minus_1_new - x_minus_1_old);
                x_minus_1_old = x_minus_1_new;
                iter = iter + 1;
            end
            
            x_minus_1 = x_minus_1_new;
            p_minus_1 = p_minus_half_new + h_over_2 * obj.hamSys.hamDx(x_minus_1_new, p_minus_half_new);
            
        end

        
    end
    
end

