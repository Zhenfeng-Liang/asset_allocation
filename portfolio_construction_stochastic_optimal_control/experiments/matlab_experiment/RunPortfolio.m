mu = [0.05,0.05];
sig = [1,1];
gamma = 10;
xT = [0.01;0.01];
pT = [0;0];
T = 1;
N = 1;
rho = [1,0;0,1];
tol = 0.001;

vlist = {'x','y'};


hamilton = Hamiltonian(vlist,@a,@b, mu, rho, sig, gamma);

[flowX, flowP] =  leapfrog(T, N, xT, pT, hamilton);
%for i = 0.1:0.1:1
 %   newtonRaphsonP(tol, hamilton, i, 0, T/N)
%end
t = 0:T/N:T;

plot(t, flowX(1,:));hold on;

% exact merton 1d lognormal

%x = xT(1,1)*exp(-mu(1)/gamma*(T-t));

%plot(t, x);



















