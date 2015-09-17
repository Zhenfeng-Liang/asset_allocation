function y = c(x, rho, sig)
    %y = rho*sig*sig*((0.05-x)^2);
    syms x;
    bv = b(x, sig);
    y = {rho(1,1)*bv{1,1}*bv{1,1}, rho(1,2)*bv{1,1}*bv{2,1};...
        rho(2,1)*bv{2,1}*bv{1,1}, rho(2,2)*bv{2,1}*bv{2,1}};
end