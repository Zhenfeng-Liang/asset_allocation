classdef UtilityCalculator
% This class includes the utility system.  
% It can calculate all utility functionality given gamma and the
% utility function type 
    
    properties
        gamma;
        oneOverGamma;
        kappa;
        type;
        
        % HARA only parameters
        a;
        b;
    end
    
    methods
        function obj = UtilityCalculator(gamma, type, a, b)
        % gamma: float, coefficient for utility function
        % type: string, utility function type
            
            if nargin < 4
                a = 1.0;
                b = 1.0;
            end 
            
            obj.a = a;
            obj.b = b;
            
            obj.gamma = gamma;
            obj.oneOverGamma = 1.0 / gamma;
            obj.type = type;
            
            if(strcmp(obj.type, 'HARA') ...
                || strcmp(obj.type, 'CRRA'))
                obj.kappa = obj.oneOverGamma - 1.0;                
            elseif(strcmp(obj.type, 'CARA'))
                obj.kappa = -1.0;
            end            
        end
        
        
        function [res] = Au(obj, v)

            res = -obj.UDer2(v) / obj.UDer(v);
        end
        
        function [res] = U(obj, v)
            
            if(strcmp(obj.type, 'CARA'))
                res = -exp(-obj.gamma * v) / obj.gamma;
            elseif(strcmp(obj.type, 'CRRA'))
                res = v^(1 - obj.gamma) / (1 - obj.gamma);
            elseif(strcmp(obj.type, 'HARA'))
                res = obj.gamma / (1 - obj.gamma) ... 
                * (obj.a + obj.b / obj.gamma * v)^(1 - obj.gamma);
            end           
        end
        
        function [res] =UDer(obj, v)

            if(strcmp(obj.type, 'CARA'))
                res = exp(-obj.gamma * v);
            elseif(strcmp(obj.type, 'CRRA'))
                res = v^(-obj.gamma);
            elseif(strcmp(obj.type, 'HARA'))
                res = obj.b * (obj.a + obj.b / obj.gamma * v)^(-obj.gamma);
            end           
            
        end
        
        function [res] = UDer2(obj, v)
    
            if(strcmp(obj.type, 'CARA'))
                res = -obj.gamma * exp(-obj.gamma * v);
            elseif(strcmp(obj.type, 'CRRA'))
                res = -obj.gamma * v^(-obj.gamma - 1);
            elseif(strcmp(obj.type, 'HARA'))
                res = -obj.b^2 * (obj.a + obj.b / obj.gamma * ...
                                  v)^(-obj.gamma - 1);
            end           
        end
    
    end 
    
end
