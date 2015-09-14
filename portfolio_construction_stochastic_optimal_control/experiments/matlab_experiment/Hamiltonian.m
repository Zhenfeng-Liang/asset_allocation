classdef Hamiltonian
   properties
       drift;
       covariance;
       parameter;
       driftDx;
       covDx;
       covInverse;
   end
   methods
       %a and c are function handles
       function obj = Hamiltonian(a, c, mu, rho, sigma, gamma)
           obj.drift = a;
           obj.covariance = c;
           obj.parameter.drift = mu;
           obj.parameter.cor = rho;
           obj.parameter.std = sigma;
           obj.parameter.mass = gamma;
           syms x;
           obj.driftDx = matlabFunction(diff(obj.drift(x, obj.parameter.drift), x));
           obj.covDx = matlabFunction(diff(obj.covariance(x, obj.parameter.cor, obj.parameter.std), x));
           obj.covInverse = matlabFunction(inv(obj.covariance(x, obj.parameter.cor, obj.parameter.std)));
       end
       function val = covInverseDx(obj, position)
           cInv = Hamiltonian.evaluate(obj.covInverse, position);
           cDx = Hamiltonian.evaluate(obj.covDx, position);
           val = -1*cInv*cInv*cDx;
           val = double(val);
       end
       function val = hamDx(obj, x, p)
           gamma = obj.parameter.mass;
           kappa = 1/gamma-1;
           a = Hamiltonian.evaluate(obj.drift, x, obj.parameter.drift);
           aDx = Hamiltonian.evaluate(obj.driftDx, x);
           
           cDx = Hamiltonian.evaluate(obj.covDx, x);
           cInv = Hamiltonian.evaluate(obj.covInverse, x);
           cInvDx = obj.covInverseDx(x);
           
           val = 0.5/gamma*p*cDx*p + 1/gamma*p*aDx...
               + kappa*a*cInv*aDx + 0.5*kappa*a*cInvDx*a;
           val = double(val);
       end
       function val = hamDp(obj, x, p)
           gamma = obj.parameter.mass;
           a = feval(obj.drift, x, obj.parameter.drift);
           c = feval(obj.covariance, x, obj.parameter.cor, obj.parameter.std);
           val = 1/gamma * (c * p + a);
           val = double(val);
       end
       function val = hamDxp(obj, x, p)
           gamma = obj.parameter.mass;
           cDx = Hamiltonian.evaluate(obj.covDx, x);
           aDx = Hamiltonian.evaluate(obj.driftDx, x);
           
           val = 1/gamma * (cDx * p + aDx);
           val = double(val);
       end
   end
   methods(Static)
       function val = evaluate(funh, varargin) 
            switch nargin(funh);
                case 0
                    val = feval(funh);
                case 1
                    val = feval(funh, varargin{1});
                case 2
                    val = feval(funh, varargin{1}, varargin{2});
            end
       end
   end
end















