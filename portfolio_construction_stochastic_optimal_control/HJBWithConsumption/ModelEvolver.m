classdef ModelEvolver
    
    methods        
        function [res] = EvolveEuler(obj, xCurr, t, T, rebTimeStep, corr, model)
            
        % Assume p = F
            
            numSteps = ceil((T - t) / rebTimeStep);
            rebTimeStep = (T - t) / numSteps;
            
            F = length(xCurr);
            
            L = chol(corr)';
            
            Z = normrnd(0, 1, F, numSteps);
            
            % For simplification, each time step is the same.
            v = ones(1, numSteps) * sqrt(rebTimeStep);
            scale = diag(v);
            
            dZ = L * Z * scale;

            res = xCurr;            
            for i = 1:numSteps
                
                xCurr = xCurr + model.driftV(xCurr) * rebTimeStep + model.diffV(xCurr) ...
                        * dZ(:,i);
                
                xCurr = (xCurr > 0) .* xCurr;
                res = horzcat(res, xCurr);
            end
            
            
        end
        
    end    
end

