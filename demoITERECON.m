clear
clc
addpath('utils')
addpath('src')
load data/SLICE24_InVivoData
read = 448;
nspokes = 21;
nt = 71;
% nspokes = 13;
% nt = 114;
coils = 8;
kdata = permute(reshape(virtual_slice(1,:,1:nspokes*nt,:),[read,nspokes,nt,coils]),[1,2,4,3]);
imSize=224;
clear slice
clear virtual_slice
%% iter recon
FT=cell(nt,1);
for t=1:nt
    kt=k(:,(t-1)*nspokes+1:t*nspokes);
    FT{t} = dnufft(kt,imSize,'g');
end
figure,imshow3(abs(squeeze(sens)),[0,1])
%% pregenerate mask
load mask_invivo
%% LLR & GLR
b=14;
sigma =std(reshape(kdata([1:2,end-1:end],:,:,:),[],1))*3;
opts.FT = FT;
opts.maxiter =15;
opts.block_dims=[b,b];
opts.lbd = sigma;%*(b+sqrt(nt));
opts.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
opts.lowrank = 'sliding';%or 'distinct', 'sliding'

opts.mask = rmask;
tic
[recon_LLR,lout] = CS_LLR(opts,kdata,squeeze(sens));
toc

figure,imshow(abs([recon_LLR(:,:,1),recon_LLR(:,:,25),recon_LLR(:,:,50),recon_LLR(:,:,71)]),[])

%% TempTV
sigma =std(reshape(kdata([1:2,end-1:end],:,:,:),[],1));
params.FT = FT;
params.TV = TV_Temp;
params.maxiter =15;
params.lbd = 5e-3;
params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;

tic
[recon_TV,tout] = CS_tv(params,kdata,squeeze(sens));
toc
figure,imshow(abs([recon_TV(:,:,1),recon_TV(:,:,25),recon_TV(:,:,50),recon_TV(:,:,71)]),[])
%%
LTres = 13*35*3.52e-3;
HTres = 13*35*3.52e-3;
TR = 3.52e-3;
FA = 15*(pi/180);
x_aif = double(sig_tv(1:64,1));
x_aif(1:5)=mean(x_aif(1:5));
x_cor = double(sig_tv(1:64,3));
% x_cor(1:7) = mean(x_cor(1:7));
Ktrans_init=0.25;
plotflag=1;

figure
pars = ToftsModFit(x_aif,x_cor,HTres,LTres,TR,FA,Ktrans_init,plotflag);