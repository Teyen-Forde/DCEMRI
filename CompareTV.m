%% collide
%(1)CS_training, (2)spatTV, %(3)tempTV%, (3)3DTV, (4)tTV, (5sTV, (6)gTV

clear;
addpath('data')
load SLICE27_optimized
out  = cell(6,3);%noisefree, 30dB, 20dB
%% define NUFFT_GPU operator
nt=76;
FT=cell(nt,1);
imSize=288;
N=[imSize,imSize];
nspokes=14;
coils=24;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
%     wt=abs(kt)./max(abs(kt(:)));
%         FT{tt}=NUFFT(kt,1,1,0,N,2);
    FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,ones(size(kt)));
end

%% adding noise
DATA{1} = cell(1,3);
sigma = cell(1,3);
basis = cell(1,3);

DATA{1} = kdata;
sigma{1}=std(reshape(kdata([1:2,end-1:end],:,:),[],1));%/sqrt(2);

reqSNR = 30;
[noise,sigma{2}]=awgnc(kdata,reqSNR);
DATA{2} = kdata + noise;

reqSNR = 20;
[noise,sigma{3}]=awgnc(kdata,reqSNR);
DATA{3} = kdata + noise;

load outCS
M = reshape(outCS.CS,[],nt);
[~,~,V] = svd(M,0);
basis{1}=V(:,1:4);
load outCS30
M = reshape(outCS30.CS,[],nt);
[~,~,V] = svd(M,0);
basis{2}=V(:,1:4);
load outCS20
M = reshape(outCS20.CS,[],nt);
[~,~,V] = svd(M,0);
basis{3}=V(:,1:4);
%% (1)CS_training 

clear params
params.FT = FT;
params.smaps = b1;
params.TV = TV_Temp;

params.maxiter =30;
params.lbd = 1.2e-5/3;
params.dim = 1;
for i = 1 : 3
params.basis = basis{i};
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[PCs,outCS] = CS_tv(params);
toc
CS=reshape(PCs*basis{i}',imSize,imSize,nt);
out{1,i}=rot90(CS);
end

%% (2)spatTV
clear basis
clear params
basis = eye(76);
params.FT = FT;
params.smaps = b1;
params.TV = TV_Spat;
params.basis = basis;
params.maxiter =30;figure,
subplot(221),plot(abs(V(:,1)))
hold on,plot(abs(Vt(:,1)))
subplot(222),plot(abs(V(:,2)))
hold on,plot(abs(Vt(:,2)))
subplot(223),plot(abs(V(:,3)))
hold on,plot(abs(Vt(:,3)))
subplot(224),plot(abs(V(:,4)))
hold on,plot(abs(Vt(:,4)))
params.lbd = 1.2e-5/3;
params.dim = 2;
for i = 1 : 3
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[PCCp,outCS] = CS_tv(params);
toc
CS=reshape(PCCp*basis',imSize,imSize,nt);
out{2,i}=rot90(CS);
end

%% (3)3DTV
clear params
params.FT = FT;
params.smaps = b1;
params.TV = TV_3D;
params.basis = basis;
params.maxiter =30;
params.lbd = 1.2e-5/3;
params.dim = 3;
for i = 1 : 3
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[PCCp,outCS] = CS_tv(params);
toc
CS=reshape(PCCp*basis',imSize,imSize,nt);
out{3,i}=rot90(CS);
end
%% (4)tTV
clear params
params.FT = FT;
params.smaps = b1;
params.maxiter =20;
params.mus = 0;
params.mut = 1.2e-5;
for i = 1 : 3
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[gtv,~] = CS_gtv(params);
toc
out{4,i}=rot90(gtv);
end

%% (5)sTV
clear params
params.FT = FT;
params.smaps = b1;
params.maxiter =20;
params.mus = 1.2e-5/3;
params.mut = 0;
for i = 1 : 3
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[gtv,~] = CS_gtv(params);
toc
out{5,i}=rot90(gtv);
end
%% (6)gTV
clear params
params.FT = FT;
params.smaps = b1;
params.maxiter =20;
params.mus = 1.2e-5/3;
params.mut = 1.2e-5;
for i = 1 : 3
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[gtv,~] = CS_gtv(params);
toc
out{6,i}=rot90(gtv);
end
%%
addpath('kidney_phantomgeneration/')
load('mask_aifPS3corr1.mat')
mask_aif=mask_aif(:,:,27);
load('mask_cortexPS3_kid1corr1.mat')
mask_kid1=mask_cortex(:,:,27);
load('mask_cortexPS3_kid2corr1.mat')
mask_kid2=mask_cortex(:,:,27);
inds_aif=find(mask_aif);
inds_kid1=find(mask_kid1);
inds_kid2=find(mask_kid2);
%%
for j=1:3
    figure(j);hold on;
    
    for i=1:8
        phant=abs(reshape(rot90(newout{i,j},-1),[],nt));
        SI_aif=sum(phant(inds_aif,:));
        plot(SI_aif,'LineWidth',2); 
    end
    legend('GroundTruth','TempTV', 'CS+PCA','SpatTV','3DTV','TempGroup','SpatGroup','3DGroup')
    hold off
    figure(j+3);hold on;
    subplot(221); hold on;title('kid1 signal intensity')
    subplot(222); hold on;
    subplot(223); hold on;title('kid2 signal intensity')
    subplot(224); hold on;
     for i=1:8
        phant=abs(reshape(rot90(newout{i,j},-1),[],nt));
        SI_kid1=sum(phant(inds_kid1,:));
        SI_kid2=sum(phant(inds_kid2,:));
        subplot(221)
        plot(SI_kid1,'LineWidth',2);
        subplot(222)
        plot(SI_kid1,'LineWidth',2);
        subplot(223)
        plot(SI_kid2,'LineWidth',2);
        subplot(224)
        plot(SI_kid2,'LineWidth',2);
     end
     legend('GroundTruth','TempTV', 'CS+PCA','SpatTV','3DTV','TempGroup','SpatGroup','3DGroup')
     hold off;
end
  %%
  
  figure(3);title('20dB'); hold on;
  phant=abs(reshape(rot90(newout{1,3},-1),[],nt));
  SI_aif=sum(phant(inds_aif,:));
  plot(SI_aif,'LineWidth',2);
  
  load outCSV20
  phant=abs(reshape(rot90(outCSV20.CS,-1),[],nt));
  SI_aif=sum(phant(inds_aif,:));
  plot(SI_aif,'LineWidth',2);
  
  phant=abs(reshape(rot90(newout{2,3},-1),[],nt));
  SI_aif=sum(phant(inds_aif,:));
  plot(SI_aif,'LineWidth',2);
  
  phant=abs(reshape(rot90(newout{4,3},-1),[],nt));
  SI_aif=sum(phant(inds_aif,:));
  plot(SI_aif,'LineWidth',2);
  
  phant=abs(reshape(rot90(newout{5,3},-1),[],nt));
  SI_aif=sum(phant(inds_aif,:));
  plot(SI_aif,'LineWidth',2);
  
  
  
  legend('GroundTruth','kt-PCA\_exact','TempTV','SpatTV','3DTV')
  hold off
  













