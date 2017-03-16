%3D data, 1st: time points, 2nd:AIF,KID1,KID2, 3nd:Noiseless, 30dB, 20dB
load data/SI_27
%% some const
Tres = 21*35*3.52e-3;
TR = 3.52e-3;
FA = 10*(pi/180);


addpath('kidney_phantomgeneration/')
load('mask_aifPS3corr1.mat')
mask_aif=mask_aif(:,:,27);
load('mask_cortexPS3_kid1corr1.mat')
mask_kid1=mask_cortex(:,:,27);
load('mask_cortexPS3_kid2corr1.mat')
mask_kid2=mask_cortex(:,:,27);
inds_aif=find(mask_aif);
inds_kid1=find(mask_kid1);
inds_kid2=find(mask_kid2);

%%
nt=76;
Ktrans_init=0.25;
plotflag=1;
c=0;
HTres = 14*35*3.52e-3;
LTres = 14*35*3.52e-3;
for i=2
    for j=2
        c=c+1;
x_aif=SI_NUFFT(:,1,j);
x_cor=SI_TempTV(:,i,j);
figure
pars = ToftsModFit(x_aif,x_cor,HTres,LTres,TR,FA,Ktrans_init,plotflag);
out(c,:)=pars(1:2);
    end
end


