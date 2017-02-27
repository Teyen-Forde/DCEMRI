filename = '/mnt/data/abhishek/lumbar/KID_PERF_7816250/020816/meas_MID01309_FID28327_AX_DYN_ARTERIAL.dat';
image_obj = mapVBVD(filename,'doaverage');
kspace = image_obj.image();%size: read, coils, views, kz
clear image_obj
[read, coils, views, kz]=size(kspace);
kdata = POCS_partialfourier(kspace);
kz = size(kdata,4);
clear kspace
kspace = ifftc(kdata,4);
clear kdata
%%
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
FT1=NUFFT(k,1,1,0,N,2);
w=abs(k(:,:))/max(abs(k(:)));

%% choose nspokes and nt and define NUFFT operator
nspokes = 8;
nt = 186;
% nspokes = 13;
% nt=114;
%     nspokes = 21;
%     nt=71;
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
%%
params = struct();
params.FT = FT;
params.TV = TV_Temp;
params.maxiter =25;
params.lbd = 1.2e-4;

out=cell(kz,1);
for zz=1:kz
    kdata = kspace(:,:,:,zz);
    kdata = permute(kdata,[1,3,2]);
    %estimate coil sensitivity map
    imcoils = FT1'*(kdata.*repmat(w,[1,1,coils]));
    ncalib = 32; % use 24 calibration lines to compute compression
    ksize = [4,4]; % kernel size
    [smaps] = doEspiritCoilSensitivity(imcoils, ncalib, ksize);
    kdata = reshape(kdata(:,1:nspokes*nt,:),[read,nspokes,nt,coils]);
    kdata = permute(kdata,[1,2,4,3]);
    sigma =std(reshape(kdata([1:2,end-1:end],:,:,:),[],1));
    params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma*2;
    tic
    [CS] = CS_tv(params,kdata,smaps);
    toc
%     figure,imshow3(rot90(abs(CS(:,:,32:40))))
    out{zz}=rot90(CS);
end