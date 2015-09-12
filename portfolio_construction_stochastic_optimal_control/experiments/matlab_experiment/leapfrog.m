function flows = leapfrog(T, N, xT, pT, rho, mu, sig, gamma)
    flows = zeros(N+1,2);
    flows(N+1,1) = xT;
    flows(N+1,2) = pT;
    
    h = (T/N);

    for i = (N):-1:1
        pHalf = newtonRaphsonP(0.0001, @hamDx, @hamDxp, flows(i+1,1), flows(i+1,2), h, rho, mu, sig, gamma);
        flows(i,1) = newtonRaphsonX(0.0001, @hamDp, @hamDxp, flows(i+1,1), pHalf, h, rho, mu, sig, gamma);
        flows(i,2) = pHalf + 0.5*h*hamDx(flows(i,1), pHalf, rho, mu, sig, gamma);
    end
end

