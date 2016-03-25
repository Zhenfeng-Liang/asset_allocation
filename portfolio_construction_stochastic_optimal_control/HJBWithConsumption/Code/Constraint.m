classdef Constraint < handle
% This is a class for constraint which will be used in the back
% test engine class. Currently, this constraint class includes
% leverage ratio constraint, consumption constraint, draw down and
% maximum return constraints.
    
    properties

        % We will get rid off our position once total return reach
        % a maximium point, or the drawdown reach a point
        maxRet;
        drawDownThreshold;
       
        getOff; % flag, true to get rid of the position
        peakRet; % record the historical maximum return in the back
                 % test 
        
        on;  % true: constraint is on, vice versa
        
        maxLR; % maximum leverage ratio
        maxCR; % maximum consumption ratio
    end
        
    methods
        
        function obj = Constraint(maxRet, drawDownThreshold, maxLR, ...
                                  on)
            obj.maxRet = maxRet;
            obj.drawDownThreshold = drawDownThreshold;
            obj.getOff = false;
            obj.peakRet = 0;
            obj.on = on;
            obj.maxLR = maxLR;
            obj.maxCR = 0.05; % by default 5% maximum consumption ratio
        end
        
        function phi = imposeLeverage(obj, currW, phi, currX)

            % Leverage constraint
            LR = abs(phi)' * currX / currW;
            
            if LR > obj.maxLR
                
                phi = phi * (obj.maxLR / LR);
            end        
        end

        function cRate = imposeConsumptionConstraint(obj, cRate, ...
                                                     rebTimeStep, currW)
            
            % Annualized consumption rate, impose maximized
            % consumption constraint.
            cRate = min(cRate * 1.0, currW * obj.maxCR) / 1.0;             
            
        end
        
        function impose(obj, currRet)

        % Impose maximum return constraint and draw down
        % constraint. 
            
            if obj.on
                if currRet > obj.peakRet
                    obj.peakRet = currRet;
                    
                    if obj.peakRet > obj.maxRet
                        obj.getOff = true;
                    end                
                end
                
                if (currRet - obj.peakRet) < obj.drawDownThreshold
                    
                    obj.getOff = true;
                end    
                                
            end            
        end    
        
    end    
end
