classdef HamiltonianSystem
    % This class includes the hamitonian system. It can calculate all the
    % relevant hamitonian 
    %   
    
    properties
        portCalc;
        utiCalc;
        
        oneOverGamma;
        kappa;
    end
    
    methods
        function obj = HamiltonianSystem(portCalc, utiCalc)
            % portCalc: portfolioCalculator, an object from PortfolioCalculator
            % gamma: coefficient for utility function
            
            obj.portCalc = portCalc;
            obj.utiCalc = utiCalc;
            %obj.gamma = gamma;
            obj.oneOverGamma = utiCalc.oneOverGamma;
            obj.kappa = utiCalc.kappa;             
        end
        
        function [val] = ham(obj, x, p)
            % Hamilton function
            % Input: x, asset price vector, p, momentum vector
            % Output: hamiltonian value
                        
            a = obj.portCalc.model.driftV(x);
            
            term1 = 0.5 * p' * obj.portCalc.instCov(x) * p;
            term2 = obj.oneOverGamma * p' * a;
            term3 = 0.5 * obj.kappa * obj.oneOverGamma * a' * (obj.portCalc.instCov(x) \ a);
            val = term1 + term2 + term3;
            
        end

        function [val] = hamDx(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: gradient of hamilton with respect to x
            
            a = obj.portCalc.model.driftV(x);
            F = length(x);
            
            term1 = zeros(F, 1);
            term2 = zeros(F, 1);
            term3 = zeros(F, 1);
            
            for i = 1:F
                term1(i) = 0.5 * p' * obj.portCalc.instCovDer(x, i) * p; 
                term2(i) = obj.oneOverGamma * p' * obj.portCalc.model.driftDer(x, i);
                term3(i) = obj.kappa * a' * obj.oneOverGamma * ...
                    (obj.portCalc.instCov(x) \ obj.portCalc.model.driftDer(x, i) + 0.5 * obj.portCalc.invInstCovDer(x, i) * a);
            end
            
            val = term1 + term2 + term3;
        end
        
        function [val] = hamDxx(obj, x, p)
            % Input: x, asset price vector, p momentum vector
            % Output: a F x F matrix, second derivative of hamilton with respect to x_ij
                        
            a = obj.portCalc.model.driftV(x);
            
            F = length(x);
            term1 = zeros(F, F);
            term2 = zeros(F, F);
            term3 = zeros(F, F);
            
            for i = 1:F
                for j = 1:F
                    term1(i,j) = 0.5 * p' * obj.portCalc.instCovDer2(x, i, j) * p;
                    
                    term2(i,j) = obj.oneOverGamma * p' * obj.portCalc.model.driftDer2(x, i, j);
                    
                    term3(i,j) = obj.kappa * obj.oneOverGamma ...
                        * (a' / obj.portCalc.instCov(x) * obj.portCalc.model.driftDer2(x, i, j) ...
                           + obj.portCalc.model.driftDer(x, j)' / obj.portCalc.instCov(x) * obj.portCalc.model.driftDer(x, i) ...
                           + a' * obj.portCalc.invInstCovDer(x, j) * obj.portCalc.model.driftDer(x, i) ...
                           + a' * obj.portCalc.invInstCovDer(x, i) * obj.portCalc.model.driftDer(x, j) ...
                           + 0.5 * a' * obj.portCalc.invInstCovDer2(x, i, j) * a);
                    
                    val = term1 + term2 + term3;
                end
            end
        end
        
        
        function [val] = hamDp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: derivative of Hamilton with respect to momentum
            
            val = obj.portCalc.instCov(x) * p ...
                  + obj.oneOverGamma * obj.portCalc.model.driftV(x);
            
        end

        function [val] = hamDpp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: second derivative of Hamilton with respect to p, momentum vectors
            
            val = obj.portCalc.instCov(x);
        end
        
        function [val] = hamDxp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: second derivative matrix of Hamilton with respect to x and p
            % Note: ith row of the result matrix is the first derivative with respect to x_i, jth column is the first derivative with respect to p_j
            
            F = length(x);
            val = zeros(F, F);

            for i = 1:F
                val(i,:) = obj.portCalc.instCovDer(x, i) * p ...
                    + obj.oneOverGamma * obj.portCalc.model.driftDer(x, i);                
            end
        end


    end
    
end

