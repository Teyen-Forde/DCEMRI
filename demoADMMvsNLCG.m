addpath('utils')
addpath('src')
%%
tmp1 = bart('traj -G -x 512 -y 32');
traj=bart( 'scale 0.5', tmp1 );
% writecfl('bart-files/traj',traj);
ksp = bart('phantom -x256 -s8 -k -t', traj);
img = bart('nufft -it -t bart-files/traj',ksp);
kspc = bart('fft -u 7',img);
sens = squeeze(bart('ecalib -k4 -r24 -m1',kspc));
%%
pht = phantom(256);
read = 512; views=32;
k = dtraj(read,views);
imSize = 256;
FT = dnufft(k,imSize,'c','s');
coils=8;
ksp = FT*(repmat(pht,[1,1,coils]).*sens);
%% nufft recon
w=repmat(abs(k)./max(abs(k(:))),[1,1,coils]);
[noise,sigma] = awgnc(ksp, 30);
kdata = ksp+noise;
recon_nufft=sos(FT'*(kdata.*w));
figure,imshow(abs(recon_nufft))
ftitle('nufft recon')
FF = dnufft(k,imSize,'c');
%% nlcg recon
%setup nlcg parameters
params.imsteps = imSize;
params.dbout = 1;%1; % print debug info
params.slices2d = 1; % set to 0 if 3D FFT
params.numcoeff = 1;
params.smaps = reshape(sens,[imSize,imSize,1,coils]);
params.smaps_tscale = ones(imSize);%1./(sum(abs(params.smaps).^2,4)+1e-13);
params.dim = [read,1,views,coils,1];
params.FF = {FF};
params.t0 = 1e-1; % - initial linesearch parameter, eg. 0.1
params.beta = 0.5; % - linesearch scaling parameter (between zero and one), eg. 0.5
params.alpha = 0.01; % - linesearch target parameter, always 0.01 for me
params.maxlsiter = 50; % - max number of line searches, eg. 100
params.minlsthresh = 1e-11; %  - minimum allowed t value, eg. 1e-6
params.lsfailtol = 20; % - times allowed to reach minimum threshold, eg. 4
params.l1Smooth = 1e-15; % - makes l1 norm smooth, always 1e-15
params.pNorm = 1; % - the value of the sparse norm, always 1
params.initimg = 0;
params.removeOS = 1; % because imsteps == steps, set to 1
cgiter = 100; params.cgiter = cgiter; % very very important (when to stop)
params.data = reshape(kdata,[size(kdata,1),size(kdata,2),1,size(kdata,3),1]);
%%
params.TVspatWeight = 5e-3;
tic
[recon_nlcg, nout,ntime] = simple_repcom_nlcg(params,1);
toc
figure,imshow(abs([recon_nlcg,10*(pht-recon_nlcg)]),[0,1]);
ftitle('nlcg recon and 10x difference')
%%
opts.FT = {FF};
opts.TV = TV_Spat;
opts.maxiter =25;
opts.lbd = 5e-4;
opts.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*sigma;
%  opts.epsilon =1e-6;
tic
[recon_admm,aout] = CS_tv(opts,kdata,sens);
toc
figure,imshow(abs([recon_admm,10*(pht-recon_nlcg)]),[0,1])
ftitle('admm recon and 10x difference')
%%
nn=size(nout,2);
rel_nlcg=size(nn,1);
for n=1:nn
    rel_nlcg(n)=isame(pht(:),nout(:,n));
end

na=size(aout.x,2);
rel_admm=size(na,1);
for n=1:na
    rel_admm(n)=isame(pht(:),aout.x(:,n));
end
%%
figure,plot(ntime,log10(rel_nlcg),'Linewidth',3); hold on;
plot(aout.time,log10(rel_admm),'Linewidth',3);hold off;
legend('nlcg','admm')
faxis
xlabel('time/s')
ylabel('relative error in log scale')
%%
figure,plot(log10(rel_nlcg),'Linewidth',3); hold on;
plot(log10(rel_admm),'Linewidth',3);hold off;
legend('nlcg','admm')
faxis
xlabel('iter #')
ylabel('relative error in log scale')