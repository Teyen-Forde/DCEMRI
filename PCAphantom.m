load SLICE30
[nsteps,views,coils]=size(kdata);
viewsper=views/nt;
kdata=reshape(kdata,[nsteps,viewsper,nt,coils]);
kdata=permute(kdata,[1,2,4,3]);
%% train from low resolution images 110 in total
imSize=288;
nt=110;nspokes=8;
FT=cell(nt,1);
N=[imSize,imSize];
for tt=1:nt
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
    wt=abs(kt)./max(abs(kt(:)));
    FT{tt}=NUFFT(kt,wt,1,0,N,2);
end


%%
tsz=8;%trainning size
mask=[zeros(imSize-tsz,1);ones(tsz*2,1);zeros(imSize-tsz,1)];
mask=repmat(mask,[1,viewsper,coils]);
for tt=1:nt
    coilimg=FT{tt}'*(kdata(:,:,:,tt).*mask);
    img(:,:,tt)=sum(coilimg.*conj(b1),3)./sos(b1);
end

%%
M=reshape(img,[],nt);
[Ut,St,Vt]=svd(M,0);
figure,plot(diag(St))

%truncate at 3
basis=Vt(:,1:15);
%%
load basis
nt=110;
FT=cell(nt,1);
N=[288,288];
for tt=1:nt
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
    %         wt=abs(kt)./max(abs(kt(:)));
    FT{tt}=NUFFT(kt,1,1,0,N,2);
%     FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,ones(size(kt)));
end


% [nsteps,views,coils]=size(kdata);
% viewsper=views/nt;
% kdata=reshape(kdata,[nsteps,viewsper,nt,coils]);
% kdata=permute(kdata,[1,2,4,3]);

params.TV = TV_Temp;
params.FT = FT;
params.smaps = b1;
params.data =kdata;
params.basis = basis;
params.maxiter =10;
sigma=std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.lbd = 1.2e-5/1.5;


tic
[PCCp] = CS_tv(params);
toc
PCCpa=PCCp*basis';
CS=reshape(PCCpa,imSize,imSize,nt);
CS=rot90(CS);
%%
t=1:27:110;
CSd=[CS(:,:,t(1)),CS(:,:,t(2)),CS(:,:,t(3)),CS(:,:,t(4)),CS(:,:,t(5))];
figure;
imshow(abs(CSd),[])
%%
load imgt
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),CS(:,:,tt));
end

figure,plot(relerr)
title('PCA: RelErr of frames')
figure,imshow(abs([imgt(:,:,55),CS(:,:,55),5*abs(imgt(:,:,55)-CS(:,:,55))]),[])
ylabel('PCA')