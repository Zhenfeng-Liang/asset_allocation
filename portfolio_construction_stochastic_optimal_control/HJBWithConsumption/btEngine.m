classdef btEngine
    
    properties
        
        t;
        T;
        numSteps;
        rebTimeStep;
    end
    
    methods
        
        function obj = btEngine(t, T, rebTimeStep)
            
            obj.t = t;
            obj.T = T;
            obj.numSteps = ceil((T - t) / rebTimeStep);
            obj.rebTimeStep = (T - t) / obj.numSteps;
        end

        function [wVec, phiMat, cVec] = runBackTest(obj, histData, ...
                                                    wkbOptimizer, w0, turnedOnConsumtion)
        % Assuming the historical data is cleaned. Time steps match
        % the rebalance time step, for test purpose.
            
            wVec = zeros(1, obj.numSteps + 1);
            cVec = zeros(1, obj.numSteps);

            % Assuming the histData only have the stocks we need to trade
            phiMat = zeros(size(histData, 1), obj.numSteps);  
            
            wVec(1) = w0;
            
            % Probably want to change it later
            rebFlowTSRatio = 3.0;
            
            for i = 2:(obj.numSteps + 1)
                
                
                [phiMat(:,i - 1), cVec(i - 1)] = ...
                    wkbOptimizer.optimalControlStrategy(histData(:, i - 1), 0, ...
                                                        obj.rebTimeStep, ...
                                                        obj.rebTimeStep ...
                                                        / rebFlowTSRatio, 1e-6, ...
                                                        wVec(i - 1), ...
                                                        turnedOnConsumtion);
                
                wVec(i) = wVec(i - 1) - cVec(i - 1) * ...
                          obj.rebTimeStep + phiMat(:, i-1)' * ...
                          (histData(:, i) - histData(:, i-1)); 
                
            end
            
            totConsumptionMinusTotCost = sum(cVec) * obj.rebTimeStep ...
                - (wVec(1) - wVec(obj.numSteps + 1))
            
            annualizedRet = totConsumptionMinusTotCost / w0 ...
                / (obj.T - obj.t)
        end
        
    end
    
end
