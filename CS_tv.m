function [x,out] = padm_CG(opts)
%% set up parameters and operators
A       = opts.FT;
D       = opts.TV;
b       = opts.data;%dimension nx,ny,nc
smaps   = opts.smaps;
basis   = opts.basis;
imgSize = size(smaps,1);
coils   = size(smaps,3);
[ETL,numcoeff] = size(basis);
% tol     = opts.tol;
% proj_method= opts.proj_method;
maxiter = opts.maxiter;
epsilon = opts.epsilon;%depending on noise std deviation
% mu      = opts.mu;%Augmented Lagrangian penalty parameter
% gamma =1;
%normalize data
dscale  = norm(abs(b(:)));
% %  dscale = 1;
b       = b / dscale;
% x0      = x0/ dscale;
epsilon = epsilon / dscale;
if isfield(opts,'lbd')
    if length(opts.lbd)>1
        lbd(1,1,:) = opts.lbd;
        lbd = repmat(lbd,[imgSize,imgSize,1]);
    else
        lbd = opts.lbd;
    end
else
    lbd=1e-4;
end
dim  = opts.dim;
gamma = 1.618;
%%

% x=x0;
u=0;v=0;du=0;dv=0;
% Ax=0*b;
x=FT(b); Ax=F(x);
out.objerr= []; out.time = [];
converged=0;
iter=0;
tic;
while ~converged
    %% steepest decent method
    %     r1 = (u + du);
    %     r2 = (v + dv);
    %     gx  = FT(Ax-r2)+x-r1;
    %     Agx = F(gx);
    %     gxtgx = gx(:)'*gx(:);
    %
    %     stp = (gxtgx)./(Agx(:)'*Agx(:)+gxtgx);
    %     stp(isnan(stp)) = 0;
    %     x_old=x;
    % %     Ax_old=Ax;
    %     x = x - stp * gx;
    % %     Ax = Ax - stp * Agx;
    %
    %     t_old=t;
    %     t=(1+sqrt(1+4*t_old^2))/2;
    %     x=x+(t_old-1)/t*(x-x_old);
    %     Ax=F(x);
    %     Ax=Ax+(t_old-1)/t*(Ax-Ax_old);
    r1=u+du;
    r2=v+dv;
    %         beta = 0;
    r=FT(Ax-r2)+(x-r1);
    %     if ii<11
    %         Ar=F(r);
    %         rtr=r(:)'*r(:);
    %         stp = rtr./(Ar(:)'*Ar(:)+rtr);
    %         stp(isnan(stp))=0;
    %         x=x-stp*r;
    %         Ax=Ax-stp*Ar;
    %     else
    p=-r;
    rr0=r(:)'*r(:);
    rr=rr0;
    for k=1:10
        Ap=F(p);
        alpha=r(:)'*r(:)/(Ap(:)'*Ap(:)+p(:)'*p(:));
        %         alpha(isnan(alpha)) = 0;
        x=x+alpha*p;
        Ax=Ax+alpha*Ap;
        rr_old = rr;
        r=r+alpha*(FT(Ap)+p);
        rr = r(:)'*r(:);
        if(rr/rr0<5e-4)
%             fprintf('Number of cg iterations is %d. \n',k);
            fprintf('cg: Iterations = %d, best residual = %14.8e\n', k, rr/rr0);
            break;
        end
        beta = rr/rr_old;
        p=-r+beta*p;
    end
    %     end
    %%
    u  = proxTV(x-du);
    v  = proxL2(Ax-dv);
    du = du- gamma*(x-u);
    dv = dv -gamma*(Ax-v);
    
   %% stop criterion
    iter=iter+1;
    stopCriterion = norm(Ax(:)-b(:));%divide by bnorm which is 1(normalized before)
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
% out = out * dscale;
x   = x * dscale;

fprintf('Number of total iterations is %d. \n',iter);
return;



    function k=F(x)
        k=zeros(size(b));
        xbt  = reshape(x * basis',[imgSize,imgSize, ETL]);
        xbts  = repmat(xbt,[1,1,1,coils]);
        for ee = 1 : ETL
            k(:,:,:,ee) = A{ee}*(squeeze(xbts(:,:,ee,:)).*smaps);
        end
    end
    function x=FT(k)
        xb = zeros(imgSize,imgSize,ETL,coils);
        for ee = 1: ETL
            xb(:,:,ee,:) = A{ee}'*k(:,:,:,ee).*conj(smaps);
        end
        xb = sum(xb,4);
        x  = reshape(xb, imgSize^2, ETL)*basis;
    end

    function u=proxTV(g)
        tau=0.249;
        %         maxINNERiter=10;
        
        g=reshape(g,imgSize,imgSize,numcoeff);
        udual=zeros(imgSize,imgSize,numcoeff,dim);
        %         udual_old = udual;
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
        u=reshape(u,imgSize^2,numcoeff);
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