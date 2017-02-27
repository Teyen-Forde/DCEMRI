function kout=POCS_partialfourier(kin)
%the goal is to solve partial fourier along partition
% input: kdata, ro x coil x views x partitions
% output: fullfilled kspace data
[nx,nc,ny,nz]=size(kin);
kz_actual=48;
kout=zeros(nx,nc,ny,kz_actual);
parfor i=1:ny
    kspace=permute(squeeze(kin(:,:,i,:)),[1,3,2]);
    kout(:,:,i,:)=pocs_pf(kspace,kz_actual);

%% debug only
end
function res=pocs_pf(kspace,kz_actual)
%pocs partial fourier
%kspace, kx,kz,coils
%Zhiyang Fu 09/14/2016
niter=16;
% kz_actual=48;
[kx,kz,kc]=size(kspace);
ksymetric=kspace(:,kz_actual-kz+1:kz,:);
img_lowres=ifft2c(zpad(ksymetric,[kx,kz_actual,kc]));
img_ph=img_lowres./abs(img_lowres);

% kzpad=zpad(kspace,[kx,kz_actual,kc]);
kzpad=cat(2,kspace,zeros(kx,kz_actual-kz,kc));
ind_unsampled=kz+1:kz_actual;
for i=1:niter
    img=ifft2c(kzpad);
    img=abs(img).*img_ph;
    Y=fft2c(img);
%     kzpad_old=kzpad;
    kzpad(:,ind_unsampled,:)=Y(:,ind_unsampled,:);
%     relerr(kzpad,kzpad_old)
end
res=permute(kzpad,[1,3,2]);


