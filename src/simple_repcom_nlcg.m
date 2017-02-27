function [u, out,time] = simple_repcom_nlcg(params,basis)
% create sensitivity encoded Fourier transforms
params.F = @(u) applyF(u,params.dim,params.FF,basis,size(basis,2),params.smaps,params.imsteps,params.slices2d);
params.FT = @(k) applyFT(k,params.dim,params.FF,basis,size(basis,2),params.smaps,params.smaps_tscale,params.imsteps,params.slices2d);
% create TV operators
params.TVspat = @(u) applyTVspat(u,params.dim,params.imsteps); % coef -> TVspace
params.TVTspat = @(w) applyTVTspat(w); % TVspace -> coef
% initialize image
if max(abs(params.initimg(:))) == 0
    u = params.FT(params.data); % initial image
else
    u = params.initimg;
end
t0 = params.t0;
beta = params.beta;
alpha = params.alpha;
maxlsiter = params.maxlsiter;
% initialize NLCG
g0 = wGradient(u,params);
du = -g0;
tctr = 0;
convergence = zeros(3,params.cgiter);
out=[];time=[];
tic;
for k = 1:params.cgiter
    % 1) Evaluate without gradient
    f0 = objective(u, du, 0, params);
    t = t0;
    % 2) Reevaluate with gradient
    targ = f0 - alpha*t*abs(g0(:)'*du(:));
    f1  =  objective(u, du, t, params);
    % 3) Line search to determine optimal step size (t)
    lsiter = 0;
    while (f1 > targ) && (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        % Sometimes the line search fails
        if t < params.minlsthresh
            tctr = tctr + 1;
            disp('t very small. trying again. this is ok occasionally.');
            break;
        end
        f1  =  objective(u, du, t, params);
    end
    % If line search fails too often, end recon early with the error (this happens, perhaps less iterations would suffice)
    % or, there is a problem with the operators, check smaps and nufft
    if tctr > params.lsfailtol
        disp(['there is a problem: t is consistently small. exiting after ',num2str(k),' iterations']); return;
    end
    % If line search fails too badly, end recon early with the error (this should not happen)
    if lsiter == maxlsiter
        disp(['there is a problem: reached max line search (?).  exiting after ',num2str(k),' iterations']); return;
    end
    % 4) Adjust the line search parameters so that the next one won't take as long
    if lsiter > 1
        t0 = t0 * beta;
    end
    if lsiter < 1
        t0 = t0 / beta;
    end
    % Record information about convergence
    convergence(1,k) = abs(f1-f0); %#ok<AGROW>
    convergence(2,k) = norm(t*du(:)); %#ok<AGROW>
    convergence(3,k) = f1; %#ok<AGROW>
    % 5) Conjugate gradient calculation (Fletcher-Reeves)
    u = u + t*du;
    g1 = wGradient(u,params);
    bk_FR = g1(:)'*g1(:)/(g0(:)'*g0(:));
    % bk_PRP = g1(:)'*(g1(:)-g0(:))/(g0(:)'*g0(:));
%     bk_HS = g1(:)'*(g1(:)-g0(:))/(du(:)'*(g1(:)-g0(:)));
%     bk_DY = g1(:)'*g1(:)/(du(:)'*(g1(:)-g0(:)));
%     bk=max(0,min(bk_HS,bk_DY));
    bk=bk_FR;
    g0 = g1;
    du =  - g1 + bk* du;
    if params.dbout == 1
        % Print info about this iteration to the terminal
        display(['cgi   ',num2str(k),'/',num2str(params.cgiter),', t=',num2str(t,4),', lsi=',num2str(lsiter)]);
    end
%     out=[out,u(:)];
    out=[out,u(:)];
    time=[time,toc];
end
% final_pc_coefs = u;
return

function res = objective(u, du, t, params)
%% calculates the objective function
p = params.pNorm;
U = u + du*t;
obj = params.F(U) - params.data;
obj = obj(:)'*obj(:); % ||Fx-y||_2^2
if length(params.TVspatWeight) == 1
    if params.TVspatWeight ~= 0
        w = params.TVspat(U);
        TVspat = (w.*conj(w)+params.l1Smooth).^(p/2); % ||TV(x)||_1
    else
        TVspat = 0;
    end
    TVspat = sum(TVspat(:)*params.TVspatWeight);
else
    if length(params.TVspatWeight) ~= params.numcoeff; error('weights don''t match coefficients'); end;
    TVspat = 0;
    for cc = 1:params.numcoeff
        Ut = U(:,:,:,cc);
        if params.TVspatWeight(cc) ~= 0
            w = params.TVspat(Ut);
            TVspatt = (w.*conj(w)+params.l1Smooth).^(p/2); % ||TV(x)||_1
        else
            TVspatt = 0;
        end
        TVspat = TVspat + sum(TVspatt(:)*params.TVspatWeight(cc));
    end
end
res = obj + TVspat;
return

function grad = wGradient(u,params)
%% computes total gradient, based on individual gradients
gradObj = gOBJ(u,params);
if length(params.TVspatWeight) == 1
    if params.TVspatWeight ~= 0
        gradTVspat = params.TVspatWeight.*gTVspat(u,params);
    else
        gradTVspat = 0;
    end
else
    if length(params.TVspatWeight) ~= params.numcoeff; error('weights don''t match coefficients'); end;
    gradTVspat = zeros(size(u));
    for cc = 1:params.numcoeff
        if params.TVspatWeight(cc) ~= 0
            gradTVspat(:,:,:,cc) = params.TVspatWeight(cc).*gTVspat(u(:,:,:,cc),params);
        end
    end
end
grad = gradObj + gradTVspat;
return

function gradObj = gOBJ(u,params)
%% computes data consistency gradient
gradObj = 2*params.FT(params.F(u) - params.data);
return

function grad = gTVspat(u,params)
%% computes total variation gradient
p = params.pNorm;
Dx = params.TVspat(u);
G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.TVTspat(G);
return

function w = applyTVspat(u,dim,imsteps)
slices = dim(2);
if length(size(u))>3
    ntimes = size(u,length(size(u)));
else
    ntimes = 1;
end
w = zeros(imsteps,imsteps,slices,ntimes,2);
w(:,:,:,:,1) = u([2:end,end],:,:,:) - u(:,:,:,:);
w(:,:,:,:,2) = u(:,[2:end,end],:,:) - u(:,:,:,:);
return

function u = applyTVTspat(w)
% transpose of the spatial TV operation
rex = w([1,1:end-1],:,:,:,1) - w(:,:,:,:,1);
rex(1,:,:,:) = -w(1,:,:,:,1);
rex(end,:,:,:) = w(end-1,:,:,:,1);
rey = w(:,[1,1:end-1],:,:,2) - w(:,:,:,:,2);
rey(:,1,:,:) = -w(:,1,:,:,2);
rey(:,end,:,:) = w(:,end-1,:,:,2);
u = (rex + rey);
return

function k = applyF(u, dim, FF, basis, numcoeff, smaps, imsteps, slices2d)
nsteps = dim(1); nparts = dim(2); viewsper = dim(3); ncoils = dim(4); ETL = dim(5);
% 1) Transfrom from PC space to image space
im_te = reshape(reshape(u,[imsteps*imsteps*nparts,numcoeff])*basis',[imsteps,imsteps,nparts,ETL]);
% 2) Sample each TE
k = zeros(nsteps,viewsper,nparts,ncoils,ETL);
for ee = 1:ETL
    z = smaps.*repmat(im_te(:,:,:,ee),[1,1,1,ncoils]);
    if ~slices2d
        z = fftshift(fft(ifftshift(z,3),[],3),3); % fft across z
    end
    for cc = 1:ncoils
        for ss = 1:nparts
            k(:,:,ss,cc,ee) = reshape(FF{ee}*z(:,:,ss,cc),[nsteps,viewsper]);
        end
    end
end
return

function u = applyFT(k, dim, FF, basis, numcoeff, smaps, smaps_tscale, imsteps, slices2d)
% nsteps = dim(1);
nparts = dim(2); % viewsper = dim(3);
ncoils = dim(4); ETL = dim(5);
% 1) Transform from k space to image space
im_te = zeros(imsteps,imsteps,nparts,ETL);
for ee = 1:ETL
    z = zeros(imsteps,imsteps,nparts,ncoils);
    for cc = 1:ncoils
        for ss = 1:nparts
            ktmp = k(:,:,ss,cc,ee);
            z(:,:,ss,cc) = FF{ee}'*ktmp(:);
        end
    end
    if ~slices2d
        z = fftshift(ifft(ifftshift(z,3),[],3),3);
    end
    z = repmat(smaps_tscale,[1,1,1,ncoils]).*conj(smaps).*z;
    im_te(:,:,:,ee) = sum(z,4);
end
% 2) Transform from image space to PC space
u = (reshape(reshape(im_te,[imsteps*imsteps*nparts,ETL])*basis,[imsteps,imsteps,nparts,numcoeff]));
return