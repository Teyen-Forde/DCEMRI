addpath('utils')
addpath('src')
load data/slice
[read,nspokes,coils,nt] = size(kdata);
imSize=256;
%% composite recon, check antomical images, and generate masks for arota, kidneys
ksp = reshape(permute(kdata,[1,2,4,3]),[read,nspokes*nt,coils]);
FF = dnufft(k,imSize,'c','s');
w=repmat(abs(k)./max(abs(k(:))),[1,1,coils]);
imcomposite = sos(FF'*(ksp.*w));

figure,imshow(flipud(imcomposite),[])
maskout=drawROI(imcomposite,1);
maskin=drawROI(imcomposite,1);
mask= maskout-maskin;
figure,imshow(mask)
mask_kid2=logical(mask);
save data/mask_bart mask_aif -append
%% iter recon
FT=cell(nt,1);
for t=1:nt
    kt=k(:,(t-1)*nspokes+1:t*nspokes);
    FT{t} = dnufft(kt,imSize,'g');
end
figure,imshow3(abs(squeeze(sens)),[0,1])
%% LLR
b=16;
sigma =std(reshape(kdata([1:2,end-1:end],:,:,:),[],1))*3;
opts.FT = FT;
opts.maxiter =10;
opts.block_dims=[b,b];
opts.lbd = sigma*(b+sqrt(nt));
opts.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;

tic
[recon_LLR,lout] = CS_LLR(opts,kdata,squeeze(sens));
toc
figure,imshow(abs([recon_LLR(:,:,1),recon_LLR(:,:,25),recon_LLR(:,:,50),recon_LLR(:,:,75)]),[])
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
figure,imshow(abs([recon_TV(:,:,1),recon_TV(:,:,25),recon_TV(:,:,50),recon_TV(:,:,75)]),[])
%% check signal curve
load data/mask_bart
sig_llr=zeros(nt,3);
sig_tv=zeros(nt,3);
%%
tmp = reshape(abs(recon_TV),[],nt);
inds=find(mask_aif);
sig_tv(:,1)=mean(tmp(inds,:),1);
inds=find(mask_kid1);
sig_tv(:,2)=mean(tmp(inds,:),1);
inds=find(mask_kid2);
sig_tv(:,3)=mean(tmp(inds,:),1);
figure,plot(sig_tv(:,1),'LineWidth',3);hold on; plot(sig_tv(:,2),'LineWidth',3);plot(sig_tv(:,3),'LineWidth',3);
legend('aif','kid1','kid2'); faxis;
ftitle('TempTV: \lambda = 5e-6')
%%
tmp = reshape(abs(recon_LLR),[],nt);
inds=find(mask_aif);
sig_llr(:,1)=mean(tmp(inds,:),1);
inds=find(mask_kid1);
sig_llr(:,2)=mean(tmp(inds,:),1);
inds=find(mask_kid2);
sig_llr(:,3)=mean(tmp(inds,:),1);
figure,plot(sig_llr(:,1),'LineWidth',3);hold on; plot(sig_llr(:,2),'LineWidth',3);plot(sig_llr(:,3),'LineWidth',3);
legend('aif','kid1','kid2'); faxis;
ftitle('GLR: \lambda = 0.25*const')