classdef hamSysCalc
    % This class includes the hamitonian system. It can calculate all the
    % relevant hamitonian 
    %   
    
    properties
        port;
        gamma;
        oneOverGamma;
        kappa;
    end
    
    methods
        function obj = hamSysCalc(port, gamma)
            % port: portfolio, an object from portConstructor
            % gamma: coefficient for utility function
            
            obj.port = port;
            obj.gamma = gamma;
            obj.oneOverGamma = 1.0 / gamma;
            obj.kappa = obj.oneOverGamma - 1;            
        end
        
        function [val] = ham(obj, x, p)
            % Hamilton function
            % Input: x, asset price vector, p, momentum vector
            % Output: hamiltonian value
                        
            a = obj.port.model.driftV(x);
            
            term1 = 0.5 * obj.oneOverGamma * p' * obj.port.instCov(x) * p;
            term2 = obj.oneOverGamma * p' * a;
            term3 = 0.5 * obj.kappa * a' * obj.port.invInstCov(x) * a;
            val = term1 + term2 + term3;
            
        end

        function [val] = hamDx(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: gradient of hamilton with respect to x
            
            a = obj.port.model.driftV(x);
            F = length(x);
            
            term1 = zeros(F, 1);
            term2 = zeros(F, 1);
            term3 = zeros(F, 1);
            
            for i = 1:F
                term1(i) = 0.5 * obj.oneOverGamma * p' * obj.port.instCovDer(x, i) * p; 
                term2(i) = obj.oneOverGamma * p' * obj.port.model.driftDer(x, i);
                term3(i) = obj.kappa * a' * (obj.port.invInstCov(x) * obj.port.model.driftDer(x, i) + 0.5 * obj.port.invInstCovDer(x, i) * a);
            end
            
            val = term1 + term2 + term3;
        end
        
        function [val] = hamDxx(obj, x, p)
            % Input: x, asset price vector, p momentum vector
            % Output: a F x F matrix, second derivative of hamilton with respect to x_ij
                        
            a = obj.port.model.driftV(x);
            
            F = length(x);
            term1 = zeros(F, F);
            term2 = zeros(F, F);
            term3 = zeros(F, F);
            
            for i = 1:F
                for j = 1:F
                    term1(i,j) = 0.5 * obj.oneOverGamma * p' * obj.port.instCovDer2(x, i, j) * p;
                    
                    term2(i,j) = obj.oneOverGamma * p' * obj.port.model.driftDer2(x, i, j);
                    
                    term3(i,j) = obj.kappa * ( a' * obj.port.invInstCov(x) * obj.port.model.driftDer2(x, i, j) ...
                        + obj.port.model.driftDer(x, j)' * obj.port.invInstCov(x) * obj.port.model.driftDer(x, i) ...
                        + a' * obj.port.invInstCovDer(x, j) * obj.port.model.driftDer(x, i) ...
                        + a' * obj.port.invInstCovDer(x, i) * obj.port.model.driftDer(x, j) ...
                        + 0.5 * a' * obj.port.invInstCovDer2(x, i, j) * a);
                    
                    val = term1 + term2 + term3;
                end
            end
        end
        
        
        function [val] = hamDp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: derivative of Hamilton with respect to momentum
            
            val = obj.oneOverGamma ...
                * (obj.port.instCov(x) * p + obj.port.model.driftV(x));
            
        end

        function [val] = hamDpp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: second derivative of Hamilton with respect to p, momentum vectors
            
            val = obj.oneOverGamma * obj.port.instCov(x);
        end
        
        function [val] = hamDxp(obj, x, p)
            % Input: x, asset price vector, p, momentum vector
            % Output: second derivative matrix of Hamilton with respect to x and p
            % Note: ith row of the result matrix is the first derivative with respect to x_i, jth column is the first derivative with respect to p_j
            
            F = length(x);
            val = zeros(F, F);
            
            for i = 1:F
                val(i,:) = obj.oneOverGamma ...
                    * (obj.port.instCovDer(x, i) * p + obj.port.model.driftDer(x, i));
            end
        end


    end
    
end

