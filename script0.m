%% perform cs, kt-pca, lps and las(later) reconstruction on phantom dataset

addpath('../GRAPPA_radial/phantom_code')
load phantom_realistic2PS3
load sensdcePS3
[nx,ny,nz,nt]=size(phantom_pre45);
%% Forward projection
% load sensdcePS3.mat
read=nx*2;coils=24;kz=42;
% nt=27;nspokes=28;
views=8*110;
imSize=288;
N=[imSize,imSize];

disp('design radial trajectory')
k = complex(zeros(read,views),zeros(read,views));
% traj=linspace(-1/2,1/2,read);
traj=(-read/2+0.5:1:read/2-0.5)/read;
% deg_skip = pi / nviews;
deg_skip = (sqrt(5)-1)/2*pi;
theta = (0:views-1)*deg_skip;
% ctr = 1;
for i=1:views
    for j=1:read
        kx = cos( theta(i) ) * traj(j);
        ky = sin( theta(i) ) * traj(j);
        k(j,i) = complex(kx,ky);
    end
end

%%
nspokes_sub=8;nt_sub=110;
kdata=zeros(read,nspokes_sub,coils,nt_sub,kz);
nz=kz;
FT=cell(nt_sub,1);

for tt=1:nt_sub
    kt=k(:,(tt-1)*nspokes_sub+1:tt*nspokes_sub);
    wt=abs(kt)./max(abs(kt(:)));
    FT{tt}=NUFFT(kt,wt,1,0,N,2);
    for zz=1:nz
        %%forward NUFFT
        
        kdata(:,:,:,tt,zz)=FT{tt}*(squeeze(sens(:,:,zz,:)).*repmat(phantom_pre45(:,:,zz,tt),[1,1,coils]));
    end
end

%%
L=rot90(L); 
t=1:27:110;
Ld=[L(:,:,t(1)),L(:,:,t(2)),L(:,:,t(3)),L(:,:,t(4)),L(:,:,t(5))];
figure,imshow(abs(Ld),[]);ylabel('L')
%
figure,imshow3(abs(L(:,:,1:27:100)),[])
%%
figure,imshow3(abs(S(:,:,1:27:100)),[])


