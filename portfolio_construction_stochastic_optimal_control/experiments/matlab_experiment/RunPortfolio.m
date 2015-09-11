mu = 0.05;
sig = 0.19;
gamma = 0.01;
xT = 0.1;
pT = 0;
T = 2;
N = 100;

flows = leapfrog(T, N, xT, pT, mu, sig, gamma);

t = 0:T/N:T;

plot(t, flows(:,1));hold on;

% exact merton 1d lognormal
x = xT*exp(-mu/gamma*(T-t));
plot(t, x);



















