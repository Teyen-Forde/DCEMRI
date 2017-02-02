%% collide
%(1)JCS, (2)BCS

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
%% JCS
clear params
params.FT = FT;
params.smaps = b1;
params.TV = TV_3D;
params.basis = basis;
params.maxiter =30;
params.lbd = 1.2e-5/3;
params.mu = params.lbd * 288;
params.dim = 3;
for i = 2
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[JCS,outJCS] = JCS(params);
toc

end
%% BCS
clear params
params.FT = FT;
params.smaps = b1;
params.TV = TV_3D;
params.basis = basis;
params.maxiter =30;
params.lbd = 1.2e-5/3;
params.mu = 1;
params.dim = 3;
params.r = 40;
for i = 2
params.data = DATA{i};
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma{i};

tic
[BCS, U, V, outBCS] = BCS(params);
toc
% CS=reshape(PCCp*basis',imSize,imSize,nt);
% out{3,i}=rot90(CS);
end