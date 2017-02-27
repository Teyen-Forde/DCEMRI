function [x,out] = CS_LLRtv(opts)
%% set up parameters and operators
A       = opts.FT;
% D       = opts.TV;
b       = opts.data;%dimension nx,ny,nc
smaps   = opts.smaps;
imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
maxiter = opts.maxiter;
epsilon = opts.epsilon;%depending on noise std deviation
mu      = opts.mu;%Augmented Lagrangian penalty parameter
%normalize data
dscale  = norm(abs(b(:)));
b       = b / dscale;
epsilon = epsilon / dscale;
gamma = 1.618;
%%
AA = @(x) adjDx(Dx(x))+adjDy(Dy(x))+FT(F(x));
maxcgiter = 16;
tol = 5e-4;
% z=0; zx=0; zy=0;
v=0; vx=0; vy=0;
M=FT(b);
out.objerr= []; out.time = [];
converged=0;
iter=0;
tic;
while ~converged
    
    %%
    zx  = softL2x(Dx(M)-vx,mu);
    zy  = softL2y(Dy(M)-vy,mu);
    z   = proxL2(F(M)-v);
    %%
    %% update M using cgsolve
    bb=FT(z+v)+adjDx(zx+vx)+adjDy(zy+vy);
    M = cgsolve(AA, bb, tol, maxcgiter);
    
    FM = F(M);
    vx = vx- gamma*(Dx(M)-zx);
    vy = vy -gamma*(Dy(M)-zy);
    v  = v - gamma*(FM-z);
    
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
function y=softL2x(x,p)
sz = size(x,1);
nx=sqrt(sum(x.^2,1));
sf=max(nx-p,0)./nx;
y = x.*repmat(sf,[sz,1]);
% y(isnan(y))=0;
end
function y=softL2y(x,p)
sz = size(x,2);
nx=sqrt(sum(x.^2,2));
sf=max(nx-p,0)./nx;
y = x.*repmat(sf,[1,sz]);
% y(isnan(y))=0;
end

function res = Dx(x)
res = x([2:end,end],:,:) - x;
end

function res = Dy(x)
res = x(:,[2:end,end],:) - x;
end

function res = Dt(x)
res = x(:,:,[2:end,end]) - x;
end

function res = adjDx(x)
res= x([1,1:end-1],:,:) - x;
res(1,:,:) = -x(1,:,:);
res(end,:,:) = x(end-1,:,:);
end
function res = adjDy(x)
res= x(:,[1,1:end-1],:) - x;
res(:,1,:) = -x(:,1,:);
res(:,end,:) = x(:,end-1,:);
end
function res = adjDt(x)
res= x(:,:,[1,1:end-1]) - x;
res(:,:,1) = -x(:,:,1);
res(:,:,end) = x(:,:,end-1);
end

% function res = norm2x(vecs)
%     N = size(vecs,1);
%     res = zeros(N,1);
%     for n = 1 : N
%         res(n)=norm(vecs(n,:),2);
%     end
% end
%
% function res = norm2y(vecs)
%     N = size(vecs,2);
%     res = zeros(N,1);
%     for n = 1 : N
%         res(n)=norm(vecs(:,n),2);
%     end
% end

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