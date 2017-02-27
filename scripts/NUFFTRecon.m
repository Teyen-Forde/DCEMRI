
clear;
addpath('data')
load SLICE27_optimized
reqSNR = 20;
[noise,sigma{2}]=awgnc(kdata,reqSNR);
kdata = kdata + noise;
%% define NUFFT_GPU operator
nt=76;
FT=cell(nt,1);
imSize=288;
N=[imSize,imSize];
nspokes=14;
coils=24;
sf=(imSize/nspokes*pi)/2;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
    wt=abs(kt)./max(abs(kt(:)))*sf;
        FT{tt}=NUFFT(kt,wt,1,0,N,2);
%     FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,wt);
end
%% TR=1.72s
img = zeros(imSize,imSize,nt);
for tt = 1:nt
    imgcoil = FT{tt}'*kdata(:,:,:,tt);
    img(:,:,tt) = sum(conj(b1).*imgcoil,3)./sos(b1);
end

figure,imshow3(abs(img))

%% TR = 0.985s

nt=133;
FT1=cell(nt,1);
imSize=288;
N=[imSize,imSize];
nspokes=8;
coils=24;
sf=(imSize/nspokes*pi)/2;
for tt=1:nt
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
    wt=abs(kt)./max(abs(kt(:)))*sf;
        FT1{tt}=NUFFT(kt,wt,1,0,N,2);
%     FT1{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,wt);
end
kdata1 = permute(reshape(permute(kdata,[1,2,4,3]),[imSize*2,nspokes, nt,coils]),[1,2,4,3]);
img1 = zeros(imSize,imSize,nt);
for tt = 1:nt
    imgcoil = FT1{tt}'*kdata1(:,:,:,tt);
    img1(:,:,tt) = sum(conj(b1).*imgcoil,3)./sos(b1);
end

figure,imshow3(abs(img1))

%% generating signal curve
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
load newout
%%
nt=76;dt=1.72;

phant=abs(reshape(rot90(newout{1,1},-1),[],nt));
SI_aif=sum(phant(inds_aif,:));
SI_kid1=sum(phant(inds_kid1,:));
SI_kid2=sum(phant(inds_kid2,:));
figure
subplot(231), hold on;plot(0:dt:dt*(nt-1),SI_aif,'LineWidth',2); 
subplot(232), hold on;plot(0:dt:dt*(nt-1),SI_kid1,'LineWidth',2); 
subplot(233), hold on;plot(0:dt:dt*(nt-1),SI_kid2,'LineWidth',2); 
nt=133;dt=0.985;
% nt=76;dt=1.72;
phant=abs(reshape(img1,[],nt));
SI_aif=sum(phant(inds_aif,:));
SI_kid1=sum(phant(inds_kid1,:));
SI_kid2=sum(phant(inds_kid2,:));
SI_aif = meanfilter(SI_aif,5);
SI_kid1 = meanfilter(SI_kid1,5);
SI_kid2 = meanfilter(SI_kid2,5);
subplot(231), plot(0:dt:dt*(nt-1),SI_aif,'LineWidth',2); 
subplot(232), plot(0:dt:dt*(nt-1),SI_kid1,'LineWidth',2); 
subplot(233), plot(0:dt:dt*(nt-1),SI_kid2,'LineWidth',2); 

nt=76;dt=1.72;

phant=abs(reshape(rot90(newout{1,1},-1),[],nt));
SI_aif=sum(phant(inds_aif,:));
SI_kid1=sum(phant(inds_kid1,:));
SI_kid2=sum(phant(inds_kid2,:));
subplot(234), hold on;plot(0:dt:dt*(nt-1),SI_aif,'LineWidth',2); 
subplot(235), hold on;plot(0:dt:dt*(nt-1),SI_kid1,'LineWidth',2); 
subplot(236), hold on;plot(0:dt:dt*(nt-1),SI_kid2,'LineWidth',2); 
% nt=133;dt=0.985;
nt=76;dt=1.72;
phant=abs(reshape(img,[],nt));
SI_aif=sum(phant(inds_aif,:));
SI_kid1=sum(phant(inds_kid1,:));
SI_kid2=sum(phant(inds_kid2,:));
SI_aif = meanfilter(SI_aif,3);
SI_kid1 = meanfilter(SI_kid1,3);
SI_kid2 = meanfilter(SI_kid2,3);
subplot(234), plot(0:dt:dt*(nt-1),SI_aif,'LineWidth',2); 
subplot(235), plot(0:dt:dt*(nt-1),SI_kid1,'LineWidth',2); 
subplot(236), plot(0:dt:dt*(nt-1),SI_kid2,'LineWidth',2); 
%%
% SI_LSA = zeros(76,3,3);
for i=3
    phant=abs(reshape(rot90(outLSA20.L+outLSA20.S,-1),[],nt));
    SI_aif=sum(phant(inds_aif,:));
    SI_kid1=sum(phant(inds_kid1,:));
    SI_kid2=sum(phant(inds_kid2,:));
    SI_LSA(:,1,i)=SI_aif;
    SI_LSA(:,2,i)=SI_kid1;
    SI_LSA(:,3,i)=SI_kid2;
end
%% showing frame 21, 24
frame=21;
for i=3
BIG=[];
BIG=[BIG,abs([newout{2,i}(:,:,frame);5*abs(newout{1,i}(:,:,frame)-newout{2,i}(:,:,frame))])];
BIG=[BIG,abs([newout{5,i}(:,:,frame);5*abs(newout{1,i}(:,:,frame)-newout{5,i}(:,:,frame))])];

end
load outLSA20
out = outLSA20.L+outLSA20.S;
BIG=[BIG,abs([out(:,:,frame);5*abs(newout{1,i}(:,:,frame)-out(:,:,frame))])];
load outCSV20
out = outCSV30.CS;
BIG=[BIG,abs([out(:,:,frame);5*abs(newout{1,i}(:,:,frame)-out(:,:,frame))])];

figure,imshow(BIG,[]); ylabel('20dB'); colorbar









