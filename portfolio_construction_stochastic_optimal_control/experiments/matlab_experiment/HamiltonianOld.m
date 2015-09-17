classdef HamiltonianOld
   properties
       drift;
       diffusion;
       parameter;
       driftDx;
       diffusionDx;
   end
   methods
       %a and c are function handles
       function obj = HamiltonianOld(a, b, mu, rho, sigma, gamma)
           obj.parameter.drift = mu;
           obj.parameter.cor = rho;
           obj.parameter.std = sigma;
           obj.parameter.mass = gamma;
           
           [nrow, ~] = size(mu);
           [~, ncol] = size(rho);
           obj.drift = cell(nrow, 1);
           obj.diffusion = cell(nrow, ncol);
           obj.driftDx = cell(nrow, 1);
           obj.diffusionDx = cell(nrow, ncol);
           
           syms x;
           driftFun = feval(a, x, obj.parameter.drift);
           driftDiff = diff(a(x, obj.parameter.drift), x);
           diffusionFun = feval(b, x, obj.parameter.std);
           diffusionDiff = diff(b(x, obj.parameter.std), x);
           
           for i=1:nrow
               obj.drift{i, 1} = matlabFunction(driftFun{i, 1});
               obj.driftDx{i, 1} = matlabFunction(driftDiff(i, 1));
               for j=1:ncol
                   obj.diffusion{i, j} = matlabFunction(diffusionFun{i, j});
                   obj.diffusionDx{i, j} = matlabFunction(diffusionDiff(i, j));
               end
           end
       end
       function val = driftTerm(obj, x)
           val = Hamiltonian.evaluateFunhCells(obj.drift, x); 
       end
       function val = driftTermDx(obj, x)
           [nrow, ~] = size(obj.parameter.drift);
           val = zeros(nrow, nrow);
           for i = 1:nrow
               val(i, i) = Hamiltonian.evaluateFunh(obj.driftDx{i,1}, x);
           end
       end
       function val = covariance(obj, x)
           rho = obj.parameter.cor;
           b = Hamiltonian.evaluateFunhCells(obj.diffusion, x);
           val = b*rho*b'; 
       end
       function val = covInverse(obj, x)
           cov = obj.covariance(x);
           val = inv(cov);
       end
       function val = covDx(obj, x)
           rho = obj.parameter.cor;
           b = Hamiltonian.evaluateFunhCells(obj.diffusion, x);
           bDx = Hamiltonian.evaluateFunhCells(obj.diffusionDx, x);
           val = bDx*rho*b' + b*rho*bDx'; 
       end
       function val = covInverseDx(obj, x)
           cInv = obj.covInverse(x);
           cDx = obj.covDx(x);
           val = -1*cInv*cInv*cDx;
       end
       function val = hamDx(obj, x, p)
           gamma = obj.parameter.mass;
           kappa = 1/gamma-1;
           a = obj.driftTerm(x);
           aDx = obj.driftTermDx(x);
           
           cDx = obj.covDx(x);
           cInv = obj.covInverse(x);
           cInvDx = obj.covInverseDx(x);
           
           val = 0.5/gamma*p'*cDx*diag(p') + 1/gamma*p'*aDx...
               + kappa*a'*cInv*aDx + 0.5*kappa*a'*cInvDx*diag(a');
           val = val';
       end
       function val = hamDp(obj, x, p)
           gamma = obj.parameter.mass;
           a = obj.driftTerm(x);
           c = obj.covariance(x);
           val = 1/gamma * (c * p + a);
       end
       function val = hamDxp(obj, x, p)
           gamma = obj.parameter.mass;
           cDx = obj.covDx(x);
           aDx = obj.driftTermDx(x);
           % Here hard coded only for n dimensional Brownian motion
           ptemp = diag(p');
           val = 1/gamma * (cDx * ptemp + aDx);
       end
   end
   methods(Static)
       function val = evaluateFunh(funh, varargin)
           switch nargin(funh);
               case 0
                   val = feval(funh);
               case 1
                   val = feval(funh, varargin{1});
               case 2
                   val = feval(funh, varargin{1}, varargin{2});
           end
       end
       function val = evaluateFunhCells(funhCells, varargin)
           [nrow, ncol] = size(funhCells);
           val = zeros(nrow, ncol);
           for i = 1:nrow
               for j = 1:ncol
                   switch nargin(funhCells{i, j});
                       case 0
                           val(i, j) = feval(funhCells{i, j});
                       case 1
                           val(i, j) = feval(funhCells{i, j}, varargin{1}(i, 1));
                       case 2
                           val(i, j) = feval(funhCells{i, j}, varargin{1}(i, 1), varargin{2}(i, 1));
                   end
               end
           end
       end
   end
end





 









