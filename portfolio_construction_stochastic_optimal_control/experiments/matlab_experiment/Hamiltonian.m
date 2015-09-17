classdef Hamiltonian
    properties
        dimension;
        parameter;
        variables;
        drift;
        driftDx;
        diffusion;
        diffusionDx;
    end
    methods
        function obj = Hamiltonian(vlist, a, b, mu, rho, sigma, gamma)
            obj.parameter.drift = mu;
            obj.parameter.cor = rho;
            obj.parameter.std = sigma;
            obj.parameter.mass = gamma;
            [~,obj.dimension] = size(mu); 
            
            obj.variables = sym(vlist);
            obj.drift = a;
            obj.diffusion = b;
            obj.driftDx = obj.gradVector(obj.drift, mu);
            obj.diffusionDx = obj.gradMatrix(obj.diffusion, sigma);
        end
        function val = hamDp(obj, x, p)
            gamma = obj.parameter.mass;
            rho = obj.parameter.cor;
            mu = obj.parameter.drift;
            sig = obj.parameter.std;
            a = cell2mat(feval(obj.drift, x, mu));
            b = cell2mat(feval(obj.diffusion, x, sig));
            c = b * rho * b';
            val = 1/gamma*(p' * c + a');
        end
        function val = hamDx(obj, x, p)
            n = obj.dimension;
            gamma = obj.parameter.mass;
            kappa = 1/gamma-1;
            rho = obj.parameter.cor;
            sig = obj.parameter.std;
            a = cell2mat(feval(obj.drift, x, p));
            % careful here, in the subs(), since obj.variables has to be a
            % row vector, we use x' instead of x
            aDx = subs(obj.driftDx, obj.variables, x');
            b = cell2mat(feval(obj.diffusion, x, sig));
            bDx = subs(obj.diffusionDx, obj.variables, x');
            c = b * rho * b';
            cInv = inv(c);
            cDx = zeros(n,n,n);
            cInvDx = zeros(n,n,n);
            for i=1:n
               cDx(:,:,i) = bDx(:,:,i)*rho*b' + b*rho*bDx(:,:,i)';
               cInvDx(:,:,i) = -cInv*cInv*cDx(:,:,i);
            end
            val = zeros(1,n);
            temp = cInv*aDx;
            for i=1:n
                val(1,n) = 0.5/gamma*p'*cDx(:,:,i)*p + 1/gamma*p'*aDx(:,i)...
                    + kappa*a'*temp(:,i) + 0.5*kappa*a'*cInvDx(:,:,i)*a;
            end
        end
        function val = hamDxp(obj, x, p)
            n = obj.dimension;
            rho = obj.parameter.cor;
            gamma = obj.parameter.mass;
            sig = obj.parameter.std;
            % careful here, in the subs(), since obj.variables has to be a
            % row vector, we use x' instead of x
            aDx = subs(obj.driftDx, obj.variables, x');
            b = cell2mat(feval(obj.diffusion, x, sig));
            bDx = subs(obj.diffusionDx, obj.variables, x');
            cDx = zeros(n,n,n);
            temp = zeros(n,n);
            for i=1:n
               cDx(:,:,i) = bDx(:,:,i)*rho*b' + b*rho*bDx(:,:,i)';
               temp = temp+cDx(:,:,i)*p(i);
            end
            val = 1/gamma * (temp + aDx);
        end
        function val = gradVector(obj, funh, input)
            vlist = obj.variables;
            n = obj.dimension;
            input = reshape(input,1,n);
            for i = 1:n
                columnVector = diff(feval(funh, vlist, input),vlist(i));
                if i==1
                    val = columnVector;
                else
                    val = horzcat(val, columnVector);
                end
            end
        end
        function val = gradMatrix(obj, funh, input)
            vlist = obj.variables;
            n = obj.dimension;
            input = reshape(input, 1, n);
            for i = 1:n
                columnMatrix = diff(feval(funh, vlist, input),vlist(i));
                if i==1
                    val = columnMatrix;
                else
                    val = horzcat(val, columnMatrix);
                end
            end
            val = reshape(val, n, n, n);
        end 
    end
%     methods(Static)
%         
%     end
end



