classdef PortfolioCalculator
    % Summary: calculate asset relationship
    % Arguments: 
    %   1: model object
    %   2: correlation matrix of different Brownian Motion which drive
    %   asset price to move. Dim: p*p
    
    properties
        model;
        corrMatr;
    
    end
    
    methods
        function obj = PortfolioCalculator(model, corrMatr)
            % model: model object created by Model class
            % corrMatr: correlation matrix among different Brownian Motion
            % which drive the asset price to move
        
            obj.model = model;
            obj.corrMatr = corrMatr;
        end    
        
        function [covx] = instCov(obj, x)
            % Input: x, asset price vector, and corrMatr, correlation
            % matrix of the underlying brownian motion 
            % Output: instantaneous covariance matrix of x
            
            dofx = obj.model.diffV(x);
            covx = dofx * obj.corrMatr * dofx';
        end
        
        function [derx] = instCovDer(obj, x, i)
            % Input: x, asset price vector, derivative with respect to ith
            % asset, corrMatr, correlation matrix between underlying
            % Brownian Motion  
            % Output: derivative of covariance matrix
            
            dofx = obj.model.diffV(x);
            ddofx = obj.model.diffDer(x, i);
            derx = ddofx * obj.corrMatr * dofx' + dofx * obj.corrMatr * ddofx'; 
            
        end

        function [derx] = instCovDer2(obj, x, i, j)
            % Input: x, asset price vector, with respect to ith and jth asset
            % Output: second derivative of covariance matrix with respect
            % to ith and jth asset 
            
            dofx = obj.model.diffV(x);
            ddofxi = obj.model.diffDer(x, i);
            ddofxj = obj.model.diffDer(x, j);
            dd2ofxij = obj.model.diffDer2(x, i, j);
            
            derx = dd2ofxij * obj.corrMatr * dofx' ...
                + ddofxi * obj.corrMatr * ddofxj' ...
                + ddofxj * obj.corrMatr * ddofxi' ...
                + dofx * obj.corrMatr * obj.model.diffDer2(x, i, j)';
        end

        function [invCovx] = invInstCov(obj, x)
            % Input: x, asset price vector
            % Output: inverse matrix of the instantaneous covariance matrix

            invCovx = obj.instCov(x) \ eye(length(x));  % for test purpose
        end

        function [derx] = invInstCovDer(obj, x, i)
            % Input: x, asset price vector, with respect to ith asset,
            % correlation matrix between underlying brownian motions 
            % Output: derivative of the inverse covariance matrix with
            % respect to x_i.  i.e. (partial C_inv) / (partial x_i) 
            
            covx = obj.instCov(x);
            derx = -1 * covx \ obj.instCovDer(x, i) / covx;
            
        end

        function [derx] = invInstCovDer2(obj, x, i, j)
            % Input: x, asset price vector, with respect to ith and jth
            % assets 
            % Output: second derivative of inverse covariance matrix with
            % respect to ith and jth asset 
            
            covx = obj.instCov(x);
            covxdi = obj.instCovDer(x, i);
            covxdj = obj.instCovDer(x, j);
            derx = covx \ (covxdi / covx * covxdj + covxdj / covx * covxdi - obj.instCovDer2(x, i, j)) / covx; 
            
        end
    
    end
    
end

