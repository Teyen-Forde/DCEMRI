function [x,out] = padm_CG(opts,b,smaps)
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
gamma = 1;
%%
% x=x0;
u=0;v=0;du=0;dv=0;
% Ax=0*b;
x=FT(b); Ax=F(x);
out.objerr= []; out.time = [];out.x=[];
converged=0;
iter=0;
% AA=@(x) x+FT(F(x));
% tol=1e-3;maxcgiter=10;
tic;
while ~converged
    
    %%
    u  = proxTV(x-du);
    v  = proxL2(Ax-dv);
    du = du- gamma*(x-u);
    dv = dv -gamma*(Ax-v);
    %%
        r1=u+du;
        r2=v+dv;
        r=FT(Ax-r2)+(x-r1);
    
        p=-r;
        %     rr0=r(:)'*r(:);
        %     rr=rr0;
        rr=r(:)'*r(:);
        for k=1:10
            Ap=F(p);
            alpha=r(:)'*r(:)/(Ap(:)'*Ap(:)+p(:)'*p(:));
            %         alpha(isnan(alpha)) = 0;
            x=x+alpha*p;
            Ax=Ax+alpha*Ap;
            rr_old = rr;
            r=r+alpha*(FT(Ap)+p);
            rr = r(:)'*r(:);
            if(rr<1e-5)
                fprintf('cg: Iterations = %d, best residual = %14.8e\n', k, sqrt(rr));
                break;
            end
            beta = rr/rr_old;
            p=-r+beta*p;
        end
%     bb = FT(v+dv)+(u+du);
%     x = cgsolve(AA, bb, x, tol, maxcgiter);
%     Ax=F(x);
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
            tmp = sum(abs(Du).^2,4).^(1/2);%norm along dim 4
            snorm = sum(tmp(:))*lbd;
            %             snorm = norm(Du(:),1)*lbd;
            if(rnorm>6.18*snorm)
                lbd = lbd * 1.618;
                fprintf('lbd * 1.618 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
            elseif (snorm>6.18*rnorm)
                lbd = lbd / 1.618;
                fprintf('lbd / 1.618 = %.3e\t, rnorm = %.3e\t, snorm = %.3e\n',lbd, rnorm, snorm);
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

%     function [x, res, iter] = cgsolve(A, b, x0, tol, maxiter, verbose)
%         
%         if (nargin < 6), verbose = 1; end
%         
%         implicit = isa(A,'function_handle');
%         
%         % x = zeros(size(b));
%         x = x0;
%         r = b;
%         d = r;
%         delta = r(:)'*r(:);
%         delta0 = b(:)'*b(:);
%         numiter = 0;
%         bestx = x;
%         bestres = sqrt(delta/delta0);
%         while ((numiter < maxiter) && (delta > tol^2*delta0))
%             
%             % q = A*d
%             if (implicit), q = A(d);  else q = A*d;  end
%             
%             alpha = delta/(d(:)'*q(:));
%             x = x + alpha*d;
%             
%             if (mod(numiter+1,50) == 0)
%                 % r = b - Aux*x
%                 if (implicit), r = b-A(x);  else r = b-A*x;  end
%             else
%                 r = r - alpha*q;
%             end
%             
%             deltaold = delta;
%             delta = r(:)'*r(:);
%             beta = delta/deltaold;
%             d = r + beta*d;
%             numiter = numiter + 1;
%             if (sqrt(delta/delta0) < bestres)
%                 bestx = x;
%                 bestres = sqrt(delta/delta0);
%             end
%             
%             if ((verbose) && (mod(numiter,50)==0))
%                 disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
%                     numiter, bestres, sqrt(delta/delta0)));
%             end
%             
%         end
%         
%         if (verbose)
%             disp(sprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres));
%         end
%         x = bestx;
%         res = bestres;
%         iter = numiter;
%     end
end