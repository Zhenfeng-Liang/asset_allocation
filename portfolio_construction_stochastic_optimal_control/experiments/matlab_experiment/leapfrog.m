function [flowX, flowP] = leapfrog(T, N, xT, pT, hamilton)
    [nrow,~] = size(xT);
    flowX = zeros(nrow, N+1);
    flowP = zeros(nrow, N+1);
    flowX(:, N+1) = xT;
    flowP(:, N+1) = pT;
    
    h = (T/N);

    for i = (N):-1:1
        pHalf = newtonRaphsonP(0.0001, hamilton, flowX(:, i+1), flowP(:,i+1), h);
        flowX(:, i) = newtonRaphsonX(0.0001, hamilton, flowX(:, i+1), pHalf, h);
        flowP(:, i) = pHalf + 0.5*h*hamilton.hamDx(flowX(:, i), pHalf)';
    end
end

