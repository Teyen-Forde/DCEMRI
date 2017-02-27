load SLICE30

nt=110;
FT=cell(nt,1);
imSize=288;
N=[288,288];
nspokes=8;
for tt=1:nt
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
            wt=abs(kt)./max(abs(kt(:)));
        FT{tt}=NUFFT(kt,1,1,0,N,2);
%     FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,ones(size(kt)));
end

[nsteps,views,coils]=size(kdata);
viewsper=views/nt;
kdata=reshape(kdata,[nsteps,viewsper,nt,coils]);
kdata=permute(kdata,[1,2,4,3]);

params.TV = TV_Temp;
params.FT = FT;
params.smaps = b1;
params.data =kdata;
params.basis = eye(nt);
params.maxiter = 1;
sigma=std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.lbd = 1.2e-5/1.5;


tic
[PCCp] = CS_tv(params);
toc

CS=reshape(PCCp,imSize,imSize,nt);
% CS=rot90(CS);
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

figure,imshow(abs([imgt(:,:,25),CS(:,:,25),5*abs(imgt(:,:,25)-CS(:,:,25))]),[])

