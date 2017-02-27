function [x,out] = iter_admm(opts,b,smaps)
%% set up parameters and operators
A       = opts.FT;
D       = opts.TV;
% b       = opts.data;%dimension nx,ny,nc
% smaps   = opts.smaps;

imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
maxiter = opts.maxiter;
epsilon = opts.epsilon;%depending on noise std deviation


%normalize data
dscale  = norm(abs(b(:)));
b       = b / dscale;
epsilon = epsilon / dscale;
if isfield(opts,'lbd')
    
    lbd = opts.lbd;
else
    lbd=1e-4;
end
gamma = 1.618;
%%
% x=x0;
% u=0;v=0;
du=0;dv=0;
% Ax=0*b;
x=FT(b); 
dims = size(x);
% Ax=F(x);
out.objerr= []; out.time = [];
converged=0;
iter=0;
AA=@(a) vec(FT(F(reshape(a,dims)))+reshape(a,dims));
tic;
while ~converged
    Ax = F(x);
    u  = proxTV(x-du);
    v  = proxL2(Ax-dv);
    du = du- gamma*(x-u);
    dv = dv -gamma*(Ax-v);
    %%
    
    bb=FT(v+dv)+(u+du);
    tmp=symmlq(AA,bb(:),1e-3,10,[],[],x(:));
    x = reshape(tmp,dims);
    %% stop criterion
    iter=iter+1;
    stopCriterion = norm(Ax(:)-b(:));%divide by bnorm which is 1(normalized before)
    objerr = dscale*stopCriterion;
    out.objerr(iter) = objerr;
    out.time(iter) = toc;
    
    fprintf('Iter %d/%d: obj error = %g\n', iter,maxiter,objerr);
    
    
    if stopCriterion < epsilon
        disp('StopCriterion Achieved') ;
        converged = 1;
    end
    
    if iter >=maxiter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
%     x_show = reshape(x,imgSize,imgSize,ETL);
%     figure(100),imshow(rot90(abs(x(:,:,50))),[]); colorbar;
end
% out = out * dscale;
x   = x * dscale;

fprintf('Number of total iterations is %d. \n',iter);
return;



    function k=F(x)
        k=zeros(size(b));
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
        udual = D*g;
        for i=1:12%8~16 are reasonable
            DDtz=D*(D'*udual+g./lbd);
            udual=(udual-tau*DDtz)./(1+tau*abs(DDtz));
        end
        dif = D'*udual.*lbd;
        u=g + dif;
        if (lbd < 1e-3 && lbd > 1e-6)
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

    function v = proxL2(v)%epsilon-radius ball centered at b
        vd = v - b;
        n_vd = norm(vd(:));
        if n_vd >= epsilon;
            vd =  vd/n_vd*epsilon;
        end
        v = b + vd;
    end
end