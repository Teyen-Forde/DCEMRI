load SLICE30

nt=110;
FT=cell(nt,1);
imSize=288;
N=[288,288];
nspokes=8;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
            wt=abs(kt)./max(abs(kt(:)));
%         FT{tt}=NUFFT(kt,1,1,0,N,2);
    FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,wt);
end

[nsteps,views,coils]=size(kdata);
viewsper=views/nt;
kdata=reshape(kdata,[nsteps,viewsper,nt,coils]);
kdata=permute(kdata,[1,2,4,3]);

params.TV = TV_Temp;
params.FT = FT;
params.smaps = b1;
params.data =kdata;
% params.basis = eye(nt);
params.maxiter = 15;
sigma=std(reshape(kdata([1:8,end-7:end],:,:),[],1));
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
params.lbd = 1.2e-5/1.5;
params.mu = 1.2e-5/1.5*288;


tic
[L,S,out] = rpcacs(params);
toc

% CS=reshape(PCCp,imSize,imSize,nt);
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

%%
%% generate kspace data
% slice 27
load EhnPht
[nx,ny,nz,nt]=size(EhnPht);
read=nx*2;coils=24;kz=42;
% nt=27;nspokes=28;
views=14*76;
imSize=288;
N=[imSize,imSize];

disp('design radial trajectory')
k = complex(zeros(read,views),zeros(read,views));
% traj=linspace(-1/2,1/2,read);
traj=(-read/2+0.5:1:read/2-0.5)/read;
% deg_skip = pi / nviews;
deg_skip = (sqrt(5)-1)/2*pi;
theta = (0:views-1)*deg_skip;
% ctr = 1;
for i=1:views
    for j=1:read
        kx = cos( theta(i) ) * traj(j);
        ky = sin( theta(i) ) * traj(j);
        k(j,i) = complex(kx,ky);
    end
end
w=abs(k)./max(abs(k(:)));
%%
load sensdcePS3
%%
nspokes_sub=14;nt_sub=76;
kdata=zeros(read,nspokes_sub,coils,nt_sub,kz);
nz=kz;
FT=cell(nt_sub,1);

for tt=1:nt_sub
    kt=k(:,(tt-1)*nspokes_sub+1:tt*nspokes_sub);
%     FT{tt}=NUFFT(kt,wt,1,0,N,2);
    FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes_sub,kt,ones(size(kt)));
    for zz=1:nz
        %%forward NUFFT
        kdata(:,:,:,tt,zz)=FT{tt}*(squeeze(sens(:,:,zz,:)).*repmat(EhnPht(:,:,zz,tt),[1,1,coils]));
    end
end