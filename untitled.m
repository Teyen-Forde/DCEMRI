clear all;clc; close all;
%ISMRM_CHALLENGE reconstructs cardiac CINE images from the ISMRM 2014 
%reconstruction challenge using ADMM with TV (total variation), LLR 
%(locally low rank) constraints
%
% Regularization parameters: 
%   'opt.lambda1'    -spatial TV, 
%   'opt.lambda2'    -temporal TV, 
%   'opt.lambda2'    -LLR, 
%   'opt.B'          -LLR patch size 

% 
%   'usePar'    - flag for using parrallel imaging. default 1
%   'overlap'   - flag for performing LLR with overlapping patches. 
%                 default 1
%   'do_shift'  - flag for performing LLR with random shifts. default 0
%   'do_plot'   - flag for plotting images. default 0
%   'max_iter'  - maximum number of iterations in ADMM. default 100
%
%   
% Authors:  Yi Guo, Terrence Jao, Sajan Lingala, Xin Miao
% Date:     03/31/2014
addpath('global_utilities')
addpath('Computing_Function');

%% Raw File Names
caseno='09';
acc=50;
ref_file=sprintf('case%sD.mat',caseno);
mask_file=sprintf('ROI_subject%s.mat',caseno);
smap_file=sprintf('smaps_subject%s_test.mat',caseno);
sampling_file=sprintf('sampling_2drandom_%s.mat',num2str(acc));

%% Initalize options from parsed inputs
usePar = 1;
overlap = 1;
do_shift = 0;
do_plot = 0;
max_iter =100;

opt = struct('useParrallel', usePar, 'overlap', overlap, 'do_shift', ...
    do_shift, 'plot', do_plot, 'itermax', max_iter);

opt.FT = @(im) fft2c( im);
opt.IFT = @(im) ifft2c( im);
%% Load Raw Data

% 1.Full kspace data: kdata_fu; 
% True image:ims
ims=opt.IFT(kdata_fu);
[kx,ky,nt,nc]=size(ims);%Arrange the data dimension as kx-ky-nt-ncoils

% 2.Sampling pattern: idx
load(sampling_file);
idx=logical(idx);
U1=repmat(idx,[1 1 1 nc]); % sampling pattern U1 the same size as data
k=kdata_fu.*repmat(idx,[1 1 1 nc]);

% 3.Sensitivity Map
load(smap_file);
sMaps=smaps;
sMaps=single(sMaps);
sMaps=reshape(sMaps,[kx,ky,1,nc]);


if opt.useParrallel
    %Clear any open interactive parrallel computing session
    try
        delete(gcp);
    catch
    end
    parpool(12);            %Setup matlab parrallel pool
end


%% -------------------------Reconstruction----------------------------------
%% Recon parameters setting
opt.B = 5;                                      % Patch Size
opt.lambda1=0;                                  % Spatial(kx,ky) TV penalty
opt.lambda2=lambdas_TV(nn_TV);                  % Temporal TV penalty
opt.lambda3=lambdas_LLR(nn_LLR);				% LLR penalty
schatten_p=0.5;

opt.p1=0.05;            % dummy variable penaltie from two splitting
opt.p2=0.05;
opt.p3=0.05;
opt.p=schatten_p;  % salton norm

opt.class='single';

opt.update=3; % plot update number

opt.U=U1>0;
kU=k(opt.U);

opt.img0=sum(conj(repmat(sMaps,[1 1 opt.size(3) 1])).*opt.IFT(k),4); % use zero-filled images as initial guess

%% ADMM Recon

[imgR,opt,time]=ADMM2_LLR_TV(kU,sMaps,opt); % run ADMM recon
imgR=double(abs(imgR));

%% Metric computation: imgR and true_im are all magnitude images

true_im= squeeze(sum(abs(ims).^2,4)).^0.5; %root of sum of square to combine coils
alpha_scale = (imgR(:)'*true_im(:))./(imgR(:)'*imgR(:));
imgR=imgR*alpha_scale;

load(mask_file);
opt.mask_ROI=repmat(mask_ROI,[1 1 nt]);
nROI=sum(mask_ROI(:));
ny1=54;ny2=147;nx1=145;nx2=268;

% MSE
error=(imgR-true_im).*opt.mask_ROI;
dominator=true_im.*opt.mask_ROI;
nerror=sqrt(sum(abs(error(:)).^2)/sum(abs(dominator(:)).^2));

% SSIM
imgR_m=imgR.*opt.mask_ROI;
true_im_m=true_im.*opt.mask_ROI;
ssim=SSIM(imgR_m(ny1:ny2,nx1:nx2,:),true_im_m(ny1:ny2,nx1:nx2,:));


%HFE
hfe=hfen(imgR(ny1:ny2,nx1:nx2,:),true_im(ny1:ny2,nx1:nx2,:));

fprintf('-------------------------------------------------------------------------------\n');
fprintf('MSE:%0.5e SSIM:%0.5e HFEN:%0.5e \n', nerror,ssim,hfe);
    

if opt.useParrallel
    %Clear any open interactive parrallel computing session
    delete(gcp);
end




