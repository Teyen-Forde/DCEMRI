function [L, S, out] = rpcacs_tv(opts)
%% set up parameters and operators
A       = opts.FT;
D       = opts.TV;
b       = opts.data;%dimension nx,ny,nc
smaps   = opts.smaps;
% basis   = opts.basis;
imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
maxiter = opts.maxiter;
epsilon = opts.epsilon;%depending on noise std deviation
dscale  = norm(abs(b(:)));
% %  dscale = 1;
b       = b / dscale;
% x0      = x0/ dscale;
epsilon = epsilon / dscale;
if isfield(opts,'mu')
    if length(opts.mu)>1
        mu(1,1,:) = opts.mu;
        mu = repmat(mu,[imgSize,imgSize,1]);
    else
        mu = opts.mu;
    end
else
    mu=1e-4;
end
lbd = opts.lbd;
dim = opts.dim;
gamma = 1.618;
%%
% addpath PROPACK;
B = @(x) FT(F(x))+x;
converged=0;
maxcgiter=10;
tol = 5e-4;
% out.L=[];out.S=[];
out.time=[];out.objerr=[];
z=0; zs=0; zl=0;
v=0; vs=0; vl=0;
L=FT(b);
% L=reshape(L,[],ETL);

% rho=1/1.2;
% mu_underline=mu/1e4;
iter=0;
tic;
% sv=5;
while ~converged
    %     Lold=L;
    %% update L and S using cgsolve
    bs=FT(z+v-F(L))+(zs+vs);
    S = cgsolve(B, bs, tol, maxcgiter);
    bl=FT(z+v-F(S))+(zl+vl);
    L = cgsolve(B, bl, tol, maxcgiter);
    
    %% update zs, zl, z    zs=proxTV(S-vs);
    if ~lbd
        zs = proxTV(S-vs);
    else
        zs=S-vs;
    end
    tmp=reshape(L-vl,[],ETL);
    [U, St, V] = svd(tmp, 'econ');
    diagSt = diag(St);
    St=diag(SoftThresh(diagSt,mu));
    zl =reshape( U * St * V',imgSize,imgSize,ETL);
    
    FLS=F(L+S);
    z=proxL2(FLS-v);
    %% update vl, vs, v
    vl=vl-gamma*(L-zl);
    vs=vs-gamma*(S-zs);
    v=v-gamma*(FLS-z);
    %     mu=min(mu*rho, mu_underline);
    
    %% stop criterion
    iter=iter+1;
    stopCriterion = norm(FLS(:)-b(:));%divide by bnorm which is 1(normalized before)
    objerr = dscale*stopCriterion;
    %     out.L    = [out.L, L(:)];
    %     out.S    = [out.S, S(:)];
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
% out = out * dscale;
% x   = x * dscale;
L=L*dscale;
S=S*dscale;
% fprintf('Number of total iterations is %d. \n',ii);
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
        %         udual=D*g;
        udual = zeros(imgSize,imgSize,ETL,dim) ;
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
                fprintf('lbd * 2 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
            elseif (snorm>6.18*rnorm)
                lbd = lbd / 1.618;
                fprintf('lbd / 2 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
            end
        end
        %         u=reshape(u,imgSize^2,numcoeff);
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
% soft-thresholding function
    function y=SoftThresh(x,p)
        y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
        y(isnan(y))=0;
    end
    function y=SoftThreshNC(x,p)%NC stands for nonconvex
        p = x.^(-0.9)*p;
        y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
        y(isnan(y))=0;
    end

% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)
%
% A - Either an NxN matrix, or a function handle.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    Default = 1.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

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
end
