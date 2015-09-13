function flows = leapfrog(T, N, xT, pT, hamilton)
    flows = zeros(N+1,2);
    flows(N+1,1) = xT;
    flows(N+1,2) = pT;
    
    h = (T/N);

    for i = (N):-1:1
        pHalf = newtonRaphsonP(0.0001, hamilton, flows(i+1,1), flows(i+1,2), h);
        flows(i,1) = newtonRaphsonX(0.0001, hamilton, flows(i+1,1), pHalf, h);
        flows(i,2) = pHalf + 0.5*h*hamilton.hamDx(flows(i,1), pHalf);
    end
end

