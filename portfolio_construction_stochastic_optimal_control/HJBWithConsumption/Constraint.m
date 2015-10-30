classdef Constraint < handle
    
    properties

        % We will get rid off our position once total return reach
        % a maximium point, or the drawdown reach a point
        maxRet;
        drawDownThreshold;
       
        getOff; % flag, true to get rid of the position
        peakRet;        
    end
        
    methods
        
        function obj = Constraint(maxRet, drawDownThreshold)
            obj.maxRet = maxRet;
            obj.drawDownThreshold = drawDownThreshold;
            obj.getOff = false;
            obj.peakRet = 0;
        end
        
        function impose(obj, currRet)
            
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
