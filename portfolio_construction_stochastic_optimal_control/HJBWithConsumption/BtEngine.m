classdef BtEngine
    
    properties
        
        t;
        T;
        numSteps;
        rebTimeStep;
        constraint;
    end
    
    methods
        
        function obj = BtEngine(t, T, rebTimeStep, constraint)
            
            obj.t = t;
            obj.T = T;
            obj.numSteps = ceil((T - t) / rebTimeStep);
            obj.rebTimeStep = (T - t) / obj.numSteps;
            obj.constraint = constraint;
        end

        function [wVec, phiMat, cVec] = runBackTest(obj, histData, ...
                                                    wkbOptimizer, w0, turnedOnConsumtion)
        % Assuming the historical data is cleaned. Time steps match
        % the rebalance time step, for test purpose.
            
            wVec = zeros(1, obj.numSteps + 1);
            cVec = zeros(1, obj.numSteps);
            cumRetVec = zeros(1, obj.numSteps + 1);
            
            % Assuming the histData only have the stocks we need to trade
            phiMat = zeros(size(histData, 1), obj.numSteps);  
            
            wVec(1) = w0;
            cumRetVec(1) = 0.0;
            
            % Probably want to change it later
            % rebFlowTSRatio = 3.0;
            
            tCurr = obj.t;
            T = obj.T;
            for i = 2:(obj.numSteps + 1)
                               
                [phiMat(:,i - 1), cVec(i - 1)] = ...
                    wkbOptimizer.optimalControlStrategy(histData(:, i - 1),  ...
                                                        tCurr, T, ...
                                                        obj.rebTimeStep, ...
                                                        1e-6, wVec(i - 1), ...
                                                        turnedOnConsumtion);
                
                tCurr = tCurr + obj.rebTimeStep;
                
                wVec(i) = wVec(i - 1) - cVec(i - 1) * ...
                          obj.rebTimeStep + phiMat(:, i-1)' * ...
                          (histData(:, i) - histData(:, i-1)); 

                cumRetVec(i) = (wVec(i) - w0) / w0;
                
                obj.constraint.impose(cumRetVec(i));
                if obj.constraint.getOff
                    
                    wVec(i+1 : obj.numSteps+1) = wVec(i);
                    cumRetVec(i+1 : obj.numSteps+1) = cumRetVec(i);                    
                    break;                        
                end                                
            end
            
            totConsumptionMinusTotCost = sum(cVec) * obj.rebTimeStep ...
                - (wVec(1) - wVec(obj.numSteps + 1));
            
            annualizedRet = totConsumptionMinusTotCost / w0 ...
                / (obj.T - obj.t)
        
            
            cumRetVec;
        end
        
    end
    
end
