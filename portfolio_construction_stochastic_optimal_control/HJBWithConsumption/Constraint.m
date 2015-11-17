classdef Constraint < handle
    
    properties

        % We will get rid off our position once total return reach
        % a maximium point, or the drawdown reach a point
        maxRet;
        drawDownThreshold;
       
        getOff; % flag, true to get rid of the position
        peakRet;
        
        on;
        
        maxLR; % maximum leverage ratio
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
        end
        
        function phi = imposeLeverage(obj, currW, phi, currX)

            % Leverage constraint
            LR = abs(phi)' * currX / currW;
            
            if LR > obj.maxLR
                
                phi = phi * (obj.maxLR / LR);
            end        
        end
        
        function impose(obj, currRet)

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
