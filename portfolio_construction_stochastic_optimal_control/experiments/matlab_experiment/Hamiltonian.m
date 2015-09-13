classdef Hamiltonian
   properties
       drift;
       covariance;
       parameter;
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
       end
       function val = driftDx(obj, position)
           syms x;
           valFun(x) = diff(obj.drift(x, obj.parameter.drift), x);
           val = valFun(position);
           val = double(val);
       end
       function val = covDx(obj, position)
           syms x;
           valFun(x) = diff(obj.covariance(x, obj.parameter.cor, obj.parameter.std), x);
           val = valFun(position);
           val = double(val);
       end
       function val = covInverse(obj, position)
           val = inv(feval(obj.covariance, position, obj.parameter.cor, obj.parameter.std));
           val = double(val);
       end
       function val = covInverseDx(obj, position)
           val = -1*obj.covInverse(position)*obj.covInverse(position)*obj.covDx(position);
           val = double(val);
       end
       function val = hamDx(obj, x, p)
           gamma = obj.parameter.mass;
           kappa = 1/gamma-1;
           a = feval(obj.drift, x, obj.parameter.drift);
           aDx = obj.driftDx(x);
           cDx = obj.covDx(x);
           cInv = obj.covInverse(x);
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
           cDx = obj.covDx(x);
           aDx = obj.driftDx(x);
           
           val = 1/gamma * (cDx * p + aDx);
           val = double(val);
       end
   end
end















