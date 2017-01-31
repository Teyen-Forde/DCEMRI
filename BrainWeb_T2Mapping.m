%%
addpath('./utils')
%% make phantom
nsteps = 256*2; imsteps = nsteps/2; img_size = imsteps; N = [imsteps,imsteps];
% load
[image_t2_gtp, image_pd_gtp] = make_phantom_t2_li_brief(imsteps);
phantomsliceoi = 136;
gtt2 = image_t2_gtp(:,:,phantomsliceoi);
% mask_ncsf=~(gtt2==329);
% gtt2=mask_ncsf.*gtt2;
% mask = (gtt2>0);1
mask_gray=(gtt2==83);
mask_white=(gtt2==70);
gtpd = image_pd_gtp(:,:,phantomsliceoi);
slices = 1; nparts = length(slices);
% make smaps
coils = 8; ncoils = coils;
smaps = zeros(imsteps,imsteps,1,coils);
[xx,yy] = meshgrid(linspace(-1,1,imsteps), linspace(-1,1,imsteps));
g2 = @(x0,y0) 1/sqrt(2*pi)*exp(-((xx-x0).^2)./2 - ((yy-y0).^2)./2);
aa = linspace(0,2*pi,coils+1); aa = aa(1:end-1);
for cc = 1:coils
    smaps(:,:,1,cc) = g2(2/3*cos(aa(cc)),2/3*sin(aa(cc)));
end
%normalized smaps
smaps=smaps./repmat(sos(squeeze(smaps)),[1,1,1,coils]);
% make decay
ETL = 16; etl = ETL;
esp = 8; ESP = esp;% milliseconds
echotimes = esp*(1:1:ETL);
gtte = zeros(imsteps,imsteps,ETL);
for ee = 1:ETL
    gtte(:,:,ee) = gtpd(:,:).*exp(-echotimes(ee)./gtt2(:,:));
end
%% sample the phantom to create data
viewsper = 16; viewsperecho = viewsper;% this determines the undersampling factor (256 would be nearly full sampling)
% also, this should be a power of 2
nviews = etl*viewsper;
traj = linspace(-1/2,1/2,nsteps);
[viewtable_vec,viewtable] = getBitRevViewOrder(etl,nviews); % bit-reverse ordering of angles
deg_skip = pi / nviews;
k = zeros(nsteps,viewsper,ETL);
w = zeros(nsteps,viewsper,ETL);
kdata = zeros(nsteps,viewsper,nparts,ncoils,ETL);
FT = cell(ETL,1); % each echo has a different acquisition trajectory
%
for ee = 1:ETL
    theta = viewtable(:,ee)'*deg_skip;
    % make k-space trajectory (complex)
    for j=1:viewsperecho
        for i=1:nsteps
            k(i,j,ee) = complex(cos(theta(j))*traj(i),...
                sin(theta(j))*traj(i));
        end
    end
    % create fourier transform
    tmp = k(:,:,ee);
    FT{ee} =  NUFFT(tmp,1, 1, 0,N, 2);
    w(:,:,ee) = abs(k(:,:,ee)); % weight is sometimes useful for non-iterative reconstruction
    % create the data by sampling the underlying image
    for cc = 1:ncoils
        kdata(:,:,1,cc,ee) = FT{ee}*(smaps(:,:,cc).*gtte(:,:,ee));
    end
end

%%
% the repcom process involves principal component analysis plus NLCG
% first, PCA
numcoeff = 3; % number of pc coefficients
T2range = linspace(20,350,10000);
signal_library = zeros(ETL,length(T2range));
for ee = 1:ETL
    signal_library(ee,:) = exp(-echotimes(ee)./T2range(:));
end
[U,sigma,~] = svd(signal_library*(signal_library'));
basis = U(:,1:numcoeff);
gtt3=reshape(gtte,img_size^2,ETL)*basis;
%% set up nlcg
% first the boring parameters:
clear params
params.imsteps = imsteps;
params.dbout = 1;%1; % print debug info
params.slices2d = 1; % set to 0 if 3D FFT
params.numcoeff = numcoeff;
params.smaps = smaps;
params.smaps_tscale = ones(imsteps);%1./(sum(abs(params.smaps).^2,4)+1e-13);
params.dim = [nsteps,nparts,viewsper,ncoils,ETL];
params.FF = FT;
params.t0 = 1e-1; % - initial linesearch parameter, eg. 0.1
params.beta = 0.5; % - linesearch scaling parameter (between zero and one), eg. 0.5
params.alpha = 0.01; % - linesearch target parameter, always 0.01 for me
params.maxlsiter = 50; % - max number of line searches, eg. 100
params.minlsthresh = 1e-11; %  - minimum allowed t value, eg. 1e-6
params.lsfailtol = 20; % - times allowed to reach minimum threshold, eg. 4
params.l1Smooth = 1e-15; % - makes l1 norm smooth, always 1e-15
params.pNorm = 1; % - the value of the sparse norm, always 1
params.initimg = 0;
params.removeOS = 1; % because imsteps == steps, set to 1
cgiter = 300; params.cgiter = cgiter; % very very important (when to stop)

%%
% data=d+sigma(i)*complex(n_real,n_imag);
params.data = kdata;
params.TVspatWeight = 5e-3;
tic;
[cs_pc_map, out,time] = simple_repcom_nlcg(params,basis);
% x0 = zeros(img_size^2,numcoeff);
% [cs_pc_map]=Brain_REPCOM_padm(x0,params);
toc;

%%
RelErr_T2map=[];
% out=out1;
% out=real(out);
% for kk=[50:50:size(out,2),size(out,2)]
for kk=1:size(out,2)
    cs_pc_map=out(:,kk);
cs_imte = reshape(reshape(cs_pc_map,[imsteps*imsteps*nparts,numcoeff])*basis',[imsteps,imsteps,nparts,ETL]);
% fit t2
tmp = log(abs(cs_imte));
[m,~] = linfit(echotimes,reshape(tmp,[imsteps*imsteps*nparts,ETL]));
% m = T2mapfit(echotimes,reshape(cs_imte,[imsteps*imsteps*nparts,ETL]));
cs_t2map = abs(1./reshape(m,[imsteps,imsteps,nparts]));
T2range=329;
cs_t2map(cs_t2map > max(T2range)) = max(T2range);
cs_t2map=mask.*cs_t2map;
% out=real(out);
RelErr_T2temp=norm(cs_t2map(:)-gtt2(:))./norm(gtt2(:));
RelErr_T2map=[RelErr_T2map,RelErr_T2temp]; %#ok<AGROW>
end
M=reshape(cs_pc_map,imsteps*imsteps*nparts,numcoeff);
% RelErr_M=sqrt(sum(abs((out-repmat(gtt3(:),1,size(out,2)))).^2))./norm(gtt3(:));
RelErr_M=sqrt(sum(abs((out-repmat(gtt3(:),1,size(out,2)))).^2)/numel(gtt3));
% filename=sprintf('NLCGTrouard_%0d.mat',sigma(i)*100);
% save(filename, 'out','M', 'cs_t2map','RelErr_M','RelErr_T2map','time');
%% plot T2 map
mask = (gtt2>0);
figure,
subplot(1,2,1);
imagesc(mask.*abs(cs_t2map-gtt2)); axis image off; colormap jet; set(gca,'CLim',[0,300]); colorbar;
title([num2str((408/256*imsteps)/viewsper),'x accelerated']);
subplot(1,2,2);
imagesc(mask.*gtt2); axis image off; colormap jet; set(gca,'CLim',[0,300]); colorbar;
title('ground truth');


%% dadm
clear params
params.smaps = squeeze(smaps);
params.data = squeeze(kdata);
params.FT = FT;
params.TV = TVOP3D;
params.basis = basis;
params.maxiter =30;
params.lambda = 1e-3;
params.mu =1;
params.proj_method='TV_Chambolle';
x0 = zeros(imsteps^2,numcoeff);
tic
[cs_pc_map, out,time] = dadm(x0,params);
toc
%% padm
clear params
params.smaps = squeeze(smaps);
params.data = squeeze(kdata);
params.FT = FT;
params.TV = TVOP3D;
params.basis = basis;
params.maxiter =30;
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*1e-4;
% params.mu = 0.01;
params.proj_method='TV_Chambolle';
x0 = zeros(imsteps^2,numcoeff);
tic
[cs_pc_map, out,time] = padm(x0,params);
toc
