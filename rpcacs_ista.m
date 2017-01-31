function [L,S,out] = rpcacs_ista(opts)
%% set up parameters and operators
A       = opts.FT;
D       = opts.TV;
b       = opts.data;%dimension nx,ny,nc
smaps   = opts.smaps;
% basis   = opts.basis;
imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
% [ETL,numcoeff] = size(basis);
% tol     = opts.tol;
% proj_method= opts.proj_method;
maxiter = opts.maxiter;
% epsilon = opts.epsilon;%depending on noise std deviation
% mu and lbd are correlated by 1/imSize, i.e. lbd = mu*(1/imSize);
mu      = opts.mu;%regularization weights on low rank
lbd     = opts.lbd;%regularization weights on sparsity

%normalize data
dscale  = norm(abs(b(:)));
b       = b / dscale;

M = FT(b);
[nx,ny,nt] = size(M);
S = 0;
z = 0;
iter = 0;
converged = 0;
fprintf('\n ********** L+S _ista reconstruction **********\n')
out.objerr= []; out.time = [];
tic
while ~converged
    % low-rank update
    tmp = reshape(M-S, nx*ny, nt);
    [Ut, St, Vt] = svd(tmp,0);
    St = diag(SoftThresh(diag(St),mu));
    L = reshape(Ut*St*Vt',nx,ny,nt);
    % sparse update
    z = SoftThresh(z + D*(M-L)*0.25,lbd);
    S = M-L-D'*z;
    %data consistency
    resk = F(L+S) - b;
    M = L+S -FT(resk);
    iter = iter +1;
    objerr = dscale*norm(resk(:));
    fprintf('Iter %d/%d: obj error = %g\n', iter,maxiter,objerr);
    if(iter>=maxiter) 
        converged=1; 
    end
    out.objerr   = [out.objerr, objerr];
    out.time  = [out.time, toc];
end
L = dscale*L;
S = dscale*S;
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
% soft-thresholding function
    function y=SoftThresh(x,p)
        y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
        y(isnan(y))=0;
    end
end