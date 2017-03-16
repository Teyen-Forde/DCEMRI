function kout=POCS_partialfourier(kin, kz_actual,leftover)
%the goal is to solve partial fourier along partition
% input: kdata, ro x coil x views x partitions
% output: fullfilled kspace data
%leftover: where the slices have been leftover: beginning or end
[nx,nc,ny,nz]=size(kin);
% kz_actual=48;
kout=zeros(nx,nc,ny,kz_actual);
parfor i=1:ny
    kspace=permute(squeeze(kin(:,:,i,:)),[1,3,2]);
    kout(:,:,i,:)=pocs_pf(kspace,kz_actual,leftover);
    
    %% debug only
end

function res=pocs_pf(kspace,kz_actual,leftover)
%pocs partial fourier
%kspace, kx,kz,coils
%Zhiyang Fu 09/14/2016
niter=12;
[kx,kz,kc]=size(kspace);
ksymetric=kspace(:,kz_actual-kz+1:kz,:);
img_lowres=ifft2c(zpad(ksymetric,[kx,kz_actual,kc]));
img_ph=img_lowres./abs(img_lowres);

if lower(leftover) == 'beginning'
    kzpad=cat(2,zeros(kx,kz_actual-kz,kc),kspace);
    ind_unsampled=1:kz_actual-kz;
elseif lower(leftover) == 'end'
    kzpad=cat(2,kspace,zeros(kx,kz_actual-kz,kc));
    ind_unsampled=kz+1:kz_actual;
else
    error('choice for leftover can only be: beginning or end')
end

for iter=1:niter
    img=ifft2c(kzpad);
    img=abs(img).*img_ph;
    Y=fft2c(img);
    %     kzpad_old=kzpad;
    kzpad(:,ind_unsampled,:)=Y(:,ind_unsampled,:);
    %     relerr(kzpad,kzpad_old)
end
res=permute(kzpad,[1,3,2]);

