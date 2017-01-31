% comparing pca(with different basis)
%           rpca(with different implementations, ista or admm)
% run multiple parameter configuration for different alogrithms
% evaluate their frame error plot and ktrans estimation
clear;
addpath('data')
load SLICE27_optimized
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

%% PCA
reqSNR = 30;
% kdata0=kdata;
[noise,sigma]=awgnc(kdata,reqSNR);
data = kdata + noise;
%% trainning basis
nt=76;
FT1=cell(nt,1);
imSize=288;
N=[imSize,imSize];
nspokes=14;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
    wt=abs(kt)./max(abs(kt(:)));
%         FT{tt}=NUFFT(kt,wt,1,0,N,2);
    FT1{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,wt);
end
 tsz = 16;
 mask=[zeros(imSize-tsz,1);ones(tsz*2,1);zeros(imSize-tsz,1)];% hamming window?
 mask=repmat(mask,[1,nspokes,coils]);
for tt=1:nt
    coilimg=FT1{tt}'*(kdata(:,:,:,tt).*mask);
    img(:,:,tt)=sum(coilimg.*conj(b1),3)./sos(b1);
end
% img = squeeze(img(:,80:135,:));
M=reshape(img,[],nt);
[Ut,St,Vt]=svd(M,0);
figure,plot(diag(St))

%% 
% load EhnPht
% EhnPht = squeeze(EhnPht(:,:,27,:));
img1 = squeeze(EhnPht(:,80:135,:));
% img1 = EhnPht;
M1=reshape(EhnPht,[],nt);
[U,SS,V]=svd(M1,0);
figure,plot(diag(SS))

figure,
subplot(221),plot(abs(V(:,1)))
hold on,plot(abs(Vt(:,1)))
subplot(222),plot(abs(V(:,2)))
hold on,plot(abs(Vt(:,2)))
subplot(223),plot(abs(V(:,3)))
hold on,plot(abs(Vt(:,3)))
subplot(224),plot(abs(V(:,4)))
hold on,plot(abs(Vt(:,4)))


figure,
subplot(221),plot(abs(V(:,5)))
hold on,plot(abs(Vt(:,5)))
subplot(222),plot(abs(V(:,6)))
hold on,plot(abs(Vt(:,6)))
subplot(223),plot(abs(V(:,7)))
hold on,plot(abs(Vt(:,7)))
subplot(224),plot(abs(V(:,8)))
hold on,plot(abs(Vt(:,8)))
%%
%truncate at 3
% basis=V(:,1:4);
basis=eye(76);
clear params
params.TV = TVOP3D;
params.FT = FT;
params.smaps = b1;
params.data = kdata;
params.basis = basis;
params.maxiter =30;
sigma = std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.lbd = 1.2e-5/3;


tic
[PCCp,outCS] = CS_tv(params);
toc
PCCpa=PCCp*basis';
CS=reshape(PCCpa,imSize,imSize,nt);
CS=rot90(CS);

%%
t=1:18:76;
CSd=[CS(:,:,t(1)),CS(:,:,t(2)),CS(:,:,t(3)),CS(:,:,t(4)),CS(:,:,t(5))];
figure;
imshow(abs(CSd),[])

load imgt27
% imgt= rot90(EhnPht);
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),CS(:,:,tt));
end

figure,plot(relerr)
title('PCA: RelErr of frames')
figure,imshow(abs([imgt(:,:,19),CS(:,:,19),5*abs(imgt(:,:,19)-CS(:,:,19))]),[])
ylabel('PCA')

%% Joint TV
clear params
params.FT = FT;
params.smaps = b1;
params.data = kdata;
params.maxiter =30;
sigma = std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.mu = 1.2e-5;


tic
[jtv,outjtv] = CS_jtv(params);
toc

jtv=rot90(jtv);

%%
t=1:18:76;
jd=[jtv(:,:,t(1)),jtv(:,:,t(2)),jtv(:,:,t(3)),jtv(:,:,t(4)),jtv(:,:,t(5))];
figure;
imshow(abs(jd),[])

load imgt27
% imgt= rot90(EhnPht);
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),jtv(:,:,tt));
end

figure,plot(relerr)
title('PCA: RelErr of frames')
figure,imshow(abs([imgt(:,:,19),jtv(:,:,19),5*abs(imgt(:,:,19)-jtv(:,:,19))]),[])
ylabel('PCA')


%% Group TV
clear params
params.FT = FT;
params.smaps = b1;
params.data = kdata;
params.maxiter =30;
sigma = std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.mu = 1.2e-5/3;


tic
[gtv,outgtv] = CS_gtv(params);
toc

gtv=rot90(gtv);

%%
t=1:18:76;
gd=[gtv(:,:,t(1)),gtv(:,:,t(2)),gtv(:,:,t(3)),gtv(:,:,t(4)),gtv(:,:,t(5))];
figure;
imshow(abs(gd),[])

load imgt27
% imgt= rot90(EhnPht);
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),gtv(:,:,tt));
end

figure,plot(relerr)
title('PCA: RelErr of frames')
figure,imshow(abs([imgt(:,:,19),gtv(:,:,19),5*abs(imgt(:,:,19)-gtv(:,:,19))]),[])
ylabel('PCA')

%% RPCA_ista
clear params
params.TV = TV_Temp;
params.FT = FT1;
params.smaps = b1;
params.data =kdata;
params.maxiter = 100;
params.lbd = 1.2e-5/1.5;
params.mu = 1.2e-5/1.5*288;


profile on
tic
[L,S,outLSI] = rpcacs_ista(params);
toc
profile viewer
profile off
L=rot90(L); 
S=rot90(S); 
LplusS=L+S;

%% RPCA_admm
clear params
params.TV = TV_3D;
params.FT = FT;
params.smaps = b1;
params.data = kdata;
params.maxiter = 20;
sigma=std(reshape(kdata([1:2,end-1:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.lbd = 1.2e-5/1.5;
% params.lbd = 0;
params.mu = 1.2e-5/1.5*288;
params.dim =3;

profile on
tic
[L,S,outLSA] = rpcacs(params);
toc
profile viewer
profile off
L=rot90(L); 
S=rot90(S); 
LplusS=L+S;

%%
t=1:18:76;
LplusSd=[LplusS(:,:,t(1)),LplusS(:,:,t(2)),LplusS(:,:,t(3)),LplusS(:,:,t(4)),LplusS(:,:,t(5))];
Ld=[L(:,:,t(1)),L(:,:,t(2)),L(:,:,t(3)),L(:,:,t(4)),L(:,:,t(5))];
Sd=[S(:,:,t(1)+1)-S(:,:,t(1)),S(:,:,t(2)+1)-S(:,:,t(2)),S(:,:,t(3)+1)-S(:,:,t(3)),S(:,:,t(4)+1)-S(:,:,t(4)),S(:,:,t(5)+1)-S(:,:,t(5))];


figure;
subplot(3,1,1),imshow(abs(LplusSd),[]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[]);ylabel('S')

%%
load imgt27
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),LplusS(:,:,tt));
end

figure,plot(relerr)

figure,imshow([imgt(:,:,55),abs(LplusS(:,:,55)),5*abs(imgt(:,:,55)-LplusS(:,:,55))],[])

%%
addpath('kidney_phantomgeneration/')
load('mask_aifPS3corr1.mat')
mask_aif=mask_aif(:,:,27);
load('mask_cortexPS3_kid1corr1.mat')
mask_kid1=mask_cortex(:,:,27);
load('mask_cortexPS3_kid2corr1.mat')
mask_kid2=mask_cortex(:,:,27);

phant=abs(reshape(EhnPht,[],nt));%rot90(LplusS,-1)
inds_aif=find(mask_aif);
SI_aif=sum(phant(inds_aif,:))./numel(inds_aif);
inds_kid=find(mask_kid1);
SI_kid1=sum(phant(inds_kid,:))./numel(inds_kid);
inds_kid=find(mask_kid2);
SI_kid2=sum(phant(inds_kid,:))./numel(inds_kid);
figure,plot(SI_aif)
hold on;plot(SI_kid1)
hold on;plot(SI_kid2)
legend('aif','kid1','kid2')
