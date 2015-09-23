classdef Model
% Summary: this class defines the model dynamic and calculate the
% relative function
% Take input from the market and model coefficients, output some model
% specific value, like, drift, diff, and their first and second
% derivatives with respect to the ith asset.
    
    properties
        modelType;
        mu;
        vol;
        lambda;
        
    end
    
    methods
        function obj = Model(param)
        % modelType: string, can be 'MeanReverting' or 'LogNormal'
        % mu: mean column vector in respective models
        % vol: diffusion coeff column vector for every Brownian Motion
        % lambda: optional, for MR model, column vector, mean reverting
        % speed for each asset
            
            obj.modelType = param.modelType;
            obj.mu = param.mu;
            obj.vol = param.vol;
            if strcmp(obj.modelType, 'MeanReverting')
                obj.lambda = param.lambda;
            end
        end
        
        function [a] = driftV(obj, x)
        % Given x, asset price, output asset dynamics drift
        % Right now, I am hard coding two dimensions drift and lambda
        % Return: drift vector
            
            if strcmp(obj.modelType, 'MeanReverting')
                % mean reverting model
                a = obj.lambda .*  (obj.mu - x);  % a is the symbol on paper
            elseif strcmp(obj.modelType, 'LogNormal')
                % lognormal model
                a = obj.mu .* x;                
            end                       
        end
        
        function [der] = driftDer(obj, x, i)
        % Input: x, asset price vector, with respect to ith asset price
        % Output: drift vector derivative with respect to x_i
            
            if strcmp(obj.modelType, 'MeanReverting')
                
                % mean reverting model
                F = length(obj.lambda);
                der = zeros(F, 1);
                der(i) = -obj.lambda(i);
                
            elseif strcmp(obj.modelType, 'LogNormal')
                % lognormal model
                F = length(obj.mu);
                der = zeros(F, 1);
                der(i) = obj.mu(i);
            end            
        end
        
        function [der2] = driftDer2(obj, x, i, j)
        % Input: x, asset price vector, with respect to ith and jth asset
        % Output: second derivative of drift vector with respect to ith
        % and jth assets 
            
        % mean reverting and lognormal models
            F = length(x);
            der2 = zeros(F, 1);
            
        end
        
        function [b] = diffV(obj, x)
        % Given x, asset price, output asset dynamics diff
        % Return: diff matrix, whose dimension is n*p
        % Note: to get a dirty implementation, we let each asset price
        % only have one Brownian Motion. So b is generated as a
        % diagonal matrix. 
            
            if strcmp(obj.modelType, 'MeanReverting')
                
                % mean reverting model
                b = diag(obj.vol);        % b is the symbol in the paper
            elseif strcmp(obj.modelType, 'LogNormal')
                % lognormal model
                b = diag(obj.vol .* x);
            end
        end
        
        function [der] = diffDer(obj, x, i)
        % Input, x, asset price, ith direction
        % Output: diff derivative matrix with respect to x_i
            
            if strcmp(obj.modelType, 'MeanReverting')
                % mean reverting model
                F = length(x);
                p = size(obj.diffV(x), 2);  % p is the number of column of diff  
                der = zeros(F, p);
                
            elseif strcmp(obj.modelType, 'LogNormal')
                % lognormal model
                F = length(obj.vol);
                vecD = zeros(F, 1);
                vecD(i) = obj.vol(i);
                der = diag(vecD);
            end
        end
        
        function [der2] = diffDer2(obj, x, i, j)
        % Input, x, asset price vector, with respect to ith and jth
        % asset derivatives 
        % Output, second derivative of diff matrix with respect to ith
        % and jth asset 
            
        % mean reverting and lognormal models
            F = length(x);
            p = size(obj.diffV(x), 2);
            der2 = zeros(F, p);
        end       
    end    
end

