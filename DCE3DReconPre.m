%% a formal script including all the steps to reconstruct a 3D DCE MRI radial data.
% clear all;clc; close all;
addpath('utils')
addpath('src')
%% load data
data_dir = '/mnt/data/abhishek/lumbar/KID_PERF_7816250/020816/';
data_name = 'meas_MID01309_FID28327_AX_DYN_ARTERIAL.dat';
data_file = [data_dir,data_name];
kspace_obj = mapVBVD(data_file,'doaverage');
kspace = kspace_obj.image();
fprintf('size of kspace data: ')
disp(size(kspace))
[read,coils,views,kz] = size(kspace);
imSize = read/2;
clear kspace_obj
%% partial fourier in kz
kz_actual = 48;
leftover = 'beginning';
disp('doing POCS to address partial fourier in kz')
tic
kspace_new = POCS_partialfourier(kspace, kz_actual, leftover);
toc
disp('done.')
figure(1),subplot(121), plot(abs(squeeze(kspace(read/2,1,20,:))),'linewidth', 3)
ftitle('before POCS')
axis([0,48,0,7e-4])
xlabel('slice #')
ylabel('kspace intensity');faxis
figure(1),subplot(122), plot(1:kz_actual,abs(squeeze(kspace_new(read/2,1,20,:))),'linewidth', 3,'color','r')
ftitle('after POCS')
xlabel('slice #');
axis([0,48,0,7e-4]);faxis
delete(gcp('nocreate'))
kz = kz_actual;
%% anatomical/composite recon using nufft
k = dtraj(read,views);
FF = dnufft(k,imSize,'c','s');
w=repmat(abs(k)./max(abs(k(:))),[1,1,coils]);

% ksp = permute(cat(4,zeros(read,coils,views,5),kspace),[1,3,2,4]);
% ksp = ifftc(ksp,4);
% imcompositez = zeros(imSize,imSize,kz);
% for nz=1:kz
%     imcompositez(:,:,nz) = sos(FF'*(ksp(:,:,:,nz).*w));
% end
% 
% figure(2),
% for nz=1:kz
%     imshow(flipud(imcompositez(:,:,nz)),[0,1e-4])
%     pause(0.05)
% end

ksp = permute(kspace_new,[1,3,2,4]);
ksp = ifftc(ksp,4);
imcompositep = zeros(imSize,imSize,kz);
disp('Reconstruct composite images using nufft')
tic
for nz=1:kz
    imcompositep(:,:,nz) = sos(FF'*(ksp(:,:,:,nz).*w));
end
toc
disp('done.')

figure(3),
for nz=1:kz
    imshow(flipud(imcompositep(:,:,nz)),[0,1e-4])
    title(num2str(nz))
    pause(0.05)
end
%compare two composite images
relerr = zeros(kz,1);
for nz = 1:kz
    relerr(nz) = isame(imcompositez(:,:,nz),imcompositep(:,:,nz));
end
figure(4),plot(relerr,'linewidth',3); axis([0,50,0, 0.05])
ftitle('relative error between composite images using POCS and zero-filling')
xlabel('slice #');faxis

figure(5),imshow3(flipud(imcompositep(:,:,38:41)),[0,1e-4],[1,4])
ftitle('slice 38 to 41, apparent abnormal tissue'); faxis;

imzf = imcompositez;
impocs = imcompositep;
save imcmp imzf impocs
% writecfl('bart-files/lumbar/imcmp',imzf,impocs)
clear imcopositez imcompositep kspace

%% coil compression using bart cc
figure(100),plot(abs(squeeze(ksp(read/2,20,4,:))),'linewidth', 3)
ftitle('profile along slice: check whether having done ifft along kz or not')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% since now we are in slices rather than partitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ksp = reshape(ksp,[1,read,views,coils,kz]);
clear kspace_new
%%
ncoils = 8;
kspv = 0*ksp;
disp('bart coil compression')
tic
for nz=1:kz
    tmp = ksp(:,:,:,:,nz);
    kspv(:,:,:,:,nz) = bart('cc -r 24 -S -p 20',tmp);%%using svd method, 8 virtual coils
end
toc

%check coil images of some slices
tmp1 = bart('traj -G -x 448 -y 1491');
traj=bart( 'scale 0.5', tmp1 );
for nz = 20
    tmp = ksp(:,:,:,:,nz);
    im = bart('nufft -i -t',traj,tmp);
    tmp = kspv(:,:,:,:,nz);
    imv =bart('nufft -i -t',traj,tmp);
end
figure(6),imshow3(abs(im),[0,1e-4],[4,5])
figure(7),imshow3(abs(imv),[0,1e-4],[4,5])
%plot error image
imgh = sos(imv);
img8=sos(imv(:,:,:,1:8));
figure(8),imshow(flipud([imgh,img8,abs(imgh-img8)*10]),[0,1e-4])
ftitle('slice 20: composite image using all coils, 8 coils, and their 10xdifference')
error=zeros(coils,1);
error(1)=isame(imgh,abs(imv(:,:,1,1)));
idx = 1;
for cc=2:coils
    imc = sos(imv(:,:,:,1:cc));
    idx = idx+1;
    error(idx)=isame(imgh,imc);
end
figure(9),plot(error,'o-','linewidth',3); 
axis([1,20,0, 1])
ftitle('sos combine')
ylabel('relative error')
xlabel('# of coils');faxis

%% estimate sens maps using bart ESPIRiT
tmp1 = bart('traj -G -x 448 -y 1491');
traj=bart( 'scale 0.5', tmp1);
kspv = kspv(:,:,:,1:ncoils,:);%truncate at 8
sens = zeros(imSize,imSize,kz,ncoils);
disp('estimate 3d sens maps')
tic
for nz = 1:kz
tmp = kspv(:,:,:,:,nz);
lowres_img = bart('nufft -i -d24:24:1 -t', traj, tmp);
lowres_ksp = bart('fft -u 7', lowres_img);
% zeropad to full size
ksp_zerop = bart('resize -c 0 224 1 224', lowres_ksp);
% ESPIRiT calibration
sens(:,:,nz,:) = bart('ecalib -m1', ksp_zerop);
end
toc
% save sens sens
writecfl('bart-files/lumbar/sens',sens)
% save kspv
writecfl('bart-files/lumbar/kspv',kspv)

%% generate mask for kids and arota
sens=readcfl('bart-files/lumbar/sens');
kspv=readcfl('bart-files/lumbar/kspv');
[~,read,views,coils,kz] = size(kspv);
imSize = read/2;
%reconstruct images withconstrast frame10~64@21 spokes
inds = 250:600;
kspv = squeeze(kspv(1,:,inds,:,:));
k = dtraj(read,views);
k = k(:,inds);
FF = dnufft(k,imSize,'c','s');
w=repmat(abs(k)./max(abs(k(:))),[1,1,coils]);
imcnst = zeros(imSize,imSize,kz);
sens = permute(sens,[1,2,4,3]);
disp('Reconstruct contrast images using nufft')
tic
for nz=1:kz
    imcnst(:,:,nz) = sum(conj(sens(:,:,:,nz)).*(FF'*(kspv(:,:,:,nz).*w)),3);
end
toc
disp('done.')