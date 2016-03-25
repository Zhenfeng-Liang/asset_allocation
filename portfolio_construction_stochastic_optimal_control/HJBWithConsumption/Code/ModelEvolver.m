classdef ModelEvolver
% This is a class to simulate historical data according to the
% models you specify in the program
    
    methods        
        function [res] = EvolveEuler(obj, xCurr, t, T, rebTimeStep, corr, model)
            
        % Assume p = F
            evolveTimes = 10;
            
            numSteps = ceil((T - t) / rebTimeStep);
            numSteps = numSteps * evolveTimes;
            
            evoTimeStep = (T - t) / numSteps;
            
            F = length(xCurr);
            
            L = chol(corr)';
            
            Z = normrnd(0, 1, F, numSteps);
            
            % For simplification, each time step is the same.
            v = ones(1, numSteps) * sqrt(evoTimeStep);
            
            scale = diag(v);
            
            dZ = L * Z * scale;

            res = xCurr;            
            for i = 1:numSteps
                
                xCurr = xCurr + model.driftV(xCurr) * evoTimeStep + model.diffV(xCurr) ...
                        * dZ(:,i);
                
                xCurr = (xCurr > 0) .* xCurr;
                
                if mod(i,10) == 0
                    res = horzcat(res, xCurr);
                end
                
            end
                        
        end
        
    end    
end

