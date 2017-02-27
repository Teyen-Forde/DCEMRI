function [x,out] = CS_LLR(opts,b,smaps)
%% set up parameters and operators
A       = opts.FT;
imgSize = size(smaps,1);
coils   = size(smaps,3);
ETL     = size(b,4);
maxiter = opts.maxiter;
block_dims = opts.block_dims;
epsilon = opts.epsilon;%depending on noise std deviation
lbd = opts.lbd;

%normalize data
dscale  = norm(abs(b(:)));
b       = b / dscale;
epsilon = epsilon / dscale;
lbd     = lbd / dscale;
gamma = 1;
%%
u=0;v=0;du=0;dv=0;
x=FT(b); Ax=F(x);
out.objerr= []; out.time = [];out.x=[];
converged=0;
iter=0;
tic;
while ~converged
    u  = llr_thresh(x-du,lbd,block_dims);
%     u  = prox_LLR(x-du,lbd,block_dims);
    v  = proxL2(Ax-dv);
    du = du- gamma*(x-u);
    dv = dv -gamma*(Ax-v);
    %%
    r1=u+du;
    r2=v+dv;
    r=FT(Ax-r2)+(x-r1);
    
    p=-r;
    rr=r(:)'*r(:);
    for k=1:10
        Ap=F(p);
        alpha=r(:)'*r(:)/(Ap(:)'*Ap(:)+p(:)'*p(:));
        x=x+alpha*p;
        Ax=Ax+alpha*Ap;
        rr_old = rr;
        r=r+alpha*(FT(Ap)+p);
        rr = r(:)'*r(:);
        if(rr<1e-5)
            fprintf('cg: iter = %d, best residual = %14.8e\n', k, sqrt(rr));
            break;
        end
        beta = rr/rr_old;
        p=-r+beta*p;
    end
    %% stop criterion
    iter=iter+1;
    stopCriterion = norm(Ax(:)-b(:));%divide by bnorm which is 1(normalized before)
    objerr = dscale*stopCriterion;
    out.x = [out.x,x(:)*dscale];
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
    figure(100),imshow(rot90(abs(x(:,:,50))),[]); colorbar;
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