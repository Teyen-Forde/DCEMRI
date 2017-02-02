function [x,out] = JCS(opts)
%% set up parameters and operators
A       = opts.FT;
D       = opts.TV;
b       = opts.data;%dimension nx,ny,nc
smaps   = opts.smaps;
imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
maxiter = opts.maxiter;
epsilon = opts.epsilon;%depending on noise std deviation
lbd     = opts.lbd;%Augmented Lagrangian penalty parameter
mu      = opts.mu;%Augmented Lagrangian penalty parameter
dim     = opts.dim;
%normalize data
dscale  = norm(abs(b(:)));
b       = b / dscale;
epsilon = epsilon / dscale;
gamma = 1;
%%
AA = @(x) 2*x+FT(F(x));
maxcgiter = 16;
tol = 5e-4;
% z=0; zx=0; zy=0;
vp=0; vq=0; vr=0;
M=FT(b);
out.objerr= []; out.time = [];
converged=0;
iter=0;
tic;
while ~converged
    
    %%
    p  = SVT(M-vp);
    q  = proxTV(M-vq);
    r  = proxL2(F(M)-vr);
    %%
    %% update M using cgsolve
    bb=FT(r+vr)+(p+vp)+(q+vq);
    M = cgsolve(AA, bb, tol, maxcgiter);
    
    FM = F(M);
    vp = vp- gamma*(M-p);
    vq = vq -gamma*(M-q);
    vr = vr -gamma*(FM-r);
    
    %% stop criterion
    iter=iter+1;
    stopCriterion = norm(FM(:)-b(:));%divide by bnorm which is 1(normalized before)
    objerr = dscale*stopCriterion;
    out.objerr = [out.objerr, objerr];
    out.time = [out.time, toc];
    
    fprintf('Iter %d/%d: obj error = %g\n', iter,maxiter,objerr);
    
    
    if stopCriterion < epsilon
        disp('StopCriterion Achieved') ;
        converged = 1;
    end
    
    if iter >=maxiter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
end

x   = M * dscale;

fprintf('Number of total iterations is %d. \n',iter);
return;



    function k=F(x)
        k=zeros(size(b));
        %         xbt  = reshape(x * basis',[imgSize,imgSize, ETL]);
        xbts  = repmat(x,[1,1,1,coils]);
        for ee = 1 : ETL
            k(:,:,:,ee) = A{ee}*(squeeze(xbts(:,:,ee,:)).*smaps);
        end
    end
    function x=FT(k)
        xb = zeros(imgSize,imgSize,ETL,coils);
        for ee = 1: ETL
            xb(:,:,ee,:) = A{ee}'*k(:,:,:,ee).*conj(smaps);
        end
        x = sum(xb,4);
        %         x  = reshape(xb, imgSize^2, ETL)*basis;
    end
    function u=proxTV(g)
        tau=0.249;
%         g=reshape(g,imgSize,imgSize,numcoeff);
        udual=zeros(imgSize,imgSize,ETL,dim);
        for i=1:16%8~16 are reasonable
            DDtz=D*(D'*udual+g./lbd);
            udual=(udual-tau*DDtz)./(1+tau*abs(DDtz));
        end
        dif = D'*udual.*lbd;
        u=g + dif;
        if (lbd < 5e-4 && lbd > 5e-6)
            Du=D*u;
            rnorm = norm(dif(:),2)^2/2;
            snorm = norm(Du(:),1)*lbd;
            if(rnorm>6.18*snorm)
                lbd = lbd * 1.618;
                fprintf('lbd * 1.618 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
            elseif (snorm>6.18*rnorm)
                lbd = lbd / 1.618;
                fprintf('lbd / 1.618 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
            end
        end
    end
    function res = SVT(M)
        M = reshape(M,[],ETL);
        [U,S,V] = svd(M,0);
        diagS = diag(S);
        S = diag((diagS-mu).*(diagS>mu));%softthresholding, S are all positive
        res = reshape(U*S*V',imgSize,imgSize,ETL);
    end
    function v = proxL2(v)%epsilon-radius ball centered at b
        vd = v - b;
        n_vd = norm(vd(:));
        if n_vd >= epsilon;
            %vd = vd;
            %         else
            vd =  vd/n_vd*epsilon;
        end
        v = b + vd;
    end

end




function [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)

if (nargin < 5), verbose = 1; end

implicit = isa(A,'function_handle');

x = zeros(size(b));
r = b;
d = r;
delta = r(:)'*r(:);
delta0 = b(:)'*b(:);
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0);
while ((numiter < maxiter) & (delta > tol^2*delta0))
    
    % q = A*d
    if (implicit), q = A(d);  else q = A*d;  end
    
    alpha = delta/(d(:)'*q(:));
    x = x + alpha*d;
    
    if (mod(numiter+1,50) == 0)
        % r = b - Aux*x
        if (implicit), r = b-A(x);  else r = b-A*x;  end
    else
        r = r - alpha*q;
    end
    
    deltaold = delta;
    delta = r(:)'*r(:);
    beta = delta/deltaold;
    d = r + beta*d;
    numiter = numiter + 1;
    if (sqrt(delta/delta0) < bestres)
        bestx = x;
        bestres = sqrt(delta/delta0);
    end
    
    if ((verbose) & (mod(numiter,50)==0))
        disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
            numiter, bestres, sqrt(delta/delta0)));
    end
    
end

if (verbose)
    disp(sprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres));
end
x = bestx;
res = bestres;
iter = numiter;
end