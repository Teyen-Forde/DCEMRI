
addpath('data/')
% to reconstruct, we need kspace data, kspace trajectory and sensitivity
% map.
%%
load SLICE24_InVivoData.mat
kdata0 = permute(kdata,[1,3,2]);
kdata = kdata0;
[read, views, coils] = size(kdata);
%% generate nufft operator
k=zeros(read,views);%kspace trajectory
deg_skip = pi*((sqrt(5)-1)/2);
traj = linspace(-1/2,1/2,read);
theta=(0:views-1)*deg_skip;
for j=1:views
    for i=1:read
        k(i,j) = complex(cos(theta(j))*traj(i),sin(theta(j))*traj(i));
    end
end

imSize=224;
N=[imSize,imSize];
FT=NUFFT(k,1,1,0,N,2);
w=abs(k(:,:))/max(abs(k(:)));

%NUFFT reconstruction
imcoils = FT'*(kdata.*repmat(w,[1,1,coils]));
figure,
for cc=1:coils
    imshow(abs(imcoils(:,:,cc)),[]);
    pause(0.5)
end
% composite image
imcomposite = sqrt(sum(abs(imcoils).^2,3));% sum of square
figure,imshow(rot90(imcomposite),[])
%% smaps
ncalib = 32; % use 24 calibration lines to compute compression
ksize = [4,4]; % kernel size
[smaps] = doEspiritCoilSensitivity(imcoils, ncalib, ksize);
figure,
for cc=1:coils
    imshow(abs(smaps(:,:,cc)),[]);
    pause(0.5)
end

%% choose nspokes and nt
% nspokes = 8;
% nt = 186;
% nspokes = 13;
% nt=114;
nspokes = 21;
nt=71;
FT=cell(nt,1);
imSize=224;
N=[imSize,imSize];
% sf=(imSize/nspokes*pi/2)/2;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
%     wt=abs(kt)./max(abs(kt(:)))*sf;
%     FT{tt}=NUFFT(kt,wt,1,0,N,2);
    FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,ones(size(kt)));
end


% for tt = 1:nt
%     imgt(:,:,tt)=sum(conj(smaps).*(FT{tt}'*kdata(:,:,:,tt)),3);
% end
% figure,imshow3(rot90(abs(imgt(:,:,21:36))))
%% constants
TR=3.67e-3;
numSlices=43;
temporalresolution = TR*numSlices*nspokes
kdata = reshape(kdata(:,1:nspokes*nt,:),[read,nspokes,nt,coils]);
kdata = permute(kdata,[1,2,4,3]);
%%
params.FT = FT;
params.smaps = smaps;
params.TV = TV_3D;
params.basis = eye(nt);
params.maxiter =6;
params.lbd = 4e-6;
% params.dim = 3;

% reqSNR = 30;
% [noise,sigma]=awgnc(kdata,reqSNR);
% data = kdata + noise;

params.data = kdata;
% sigma =std(reshape(kdata([1:2,end-1:end],:,:,:),[],1))/sqrt(2);
% params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma*6;
params.epsilon = 0.078;

tic
[PCCp,outCS] = CS_tv(params);
toc
CS=reshape(PCCp,imSize,imSize,nt);
figure,imshow3(rot90(abs(CS(:,:,32:60))))

%%
clear params
params.TV = TV_Temp;
params.FT = FT;
params.smaps = smaps;
params.data = kdata;
params.maxiter = 6;
% sigma=std(reshape(kdata([1:2,end-1:end],:,:,:),[],1));
% params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma*8;
params.epsilon = 0.078;
params.lbd = 1.2e-5/1.5;
% params.lbd = 0;
params.mu = 1.2e-5/1.5*224;
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
figure,imshow3(abs(LplusS(:,:,24:39)))

%%

