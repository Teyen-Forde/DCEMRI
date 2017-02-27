% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L+S reconstruction of dynamic contrast-enhanced abdominal MRI acquired
% with golden-angle radial sampling
% 
% Temporal resolution is flexible and determined by the user in the 
% variable nspokes 
%
% Ricardo Otazo (2013)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;close all;
% addpath('nufft_toolbox/');
% addpath('utils/')
% number of spokes to be used per frame (Fibonacci number)
nspokes= 8;
% load radial data
load SLICE30.mat

[nx,ny,nc]=size(b1);
[nr,ntviews,nc]=size(kdata);
% number of frames
nt=floor(ntviews/nspokes);
% crop the data according to the number of spokes per frame
kdata=kdata(:,1:nt*nspokes,:);
k=k(:,1:nt*nspokes);
w=w(:,1:nt*nspokes);
% w = ones(size(k));
% sort the data into a time-series of undersampled images
for ii=1:nt
    kdatau(:,:,:,ii)=kdata(:,(ii-1)*nspokes+1:ii*nspokes,:);
    ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end
% multicoil NUFFT operator
param.E=MCNUFFT(ku,wu,b1);
param.d=kdatau;
recon_nufft=param.E'*param.d;
clear kdata k ku wu w
% L+S reconstruction ******************************************************
param.lambda_L=0.025;
param.lambda_S=0.5*max(abs(recon_nufft(:)));
param.nite=20;
param.tol=0.0025;
tic
[L,S] = lps_tv_ista(param);
toc
L=rot90(L); 
S=rot90(S); 
LplusS=L+S;


%%
t=1:27:110;
LplusSd=[LplusS(:,:,t(1)),LplusS(:,:,t(2)),LplusS(:,:,t(3)),LplusS(:,:,t(4)),LplusS(:,:,t(5))];
Ld=[L(:,:,t(1)),L(:,:,t(2)),L(:,:,t(3)),L(:,:,t(4)),L(:,:,t(5))];
Sd=[S(:,:,t(1)),S(:,:,t(2)),S(:,:,t(3)),S(:,:,t(4)),S(:,:,t(5))];


figure;
subplot(3,1,1),imshow(abs(LplusSd),[]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[]);ylabel('S')

Ld=Ld/pi*2;
Sd=Sd/pi*2;
LplusS=LplusS/pi*2;
%%
load imgt
relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),LplusS(:,:,tt));
end

figure,plot(relerr)

figure,imshow([imgt(:,:,55),abs(LplusS(:,:,55)),5*abs(imgt(:,:,55)-LplusS(:,:,55))],[])
%%
figure, subplot(121),imshow([imgt(:,:,55)],[]);
subplot(122),imshow(abs([LplusS(:,:,55)]),[]);




