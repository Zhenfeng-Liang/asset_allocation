mu = 0.05;
sig = 0.19;
gamma = 0.01;
xT = 0.1;
pT = 0;
T = 2;
N = 100;
rho = 1;

flows = leapfrog(T, N, xT, pT, rho, mu, sig, gamma);
%(1/gamma-1)*a(xT, mu)*cInverse(xT, rho, sig)*aDx(xT, mu)+ 0.5*(1/gamma-1)*a(xT, mu)*cInverseDx(xT, rho, sig)*a(xT,mu)
t = 0:T/N:T;

plot(t, flows(:,1));%hold on;

% exact merton 1d lognormal
%x = xT*exp(-mu/gamma*(T-t));
%plot(t, x);


















