
% This function analyzes the data using the 2-compartment curve fitting. 
% This script averages the ROI signal in the preprocessing step and 
% processes only 1 kidney.

load('phantom60.mat');

addpath('your comapartment model folder');
mask_cor = squeeze(roiMask.RCMask(:,:,:));
%mask_cor = repmat((mask_cor),[1,1,1,151]);
mask_med = squeeze(roiMask.RMMask(:,:,:));
mask_CS = squeeze(roiMask.RCSMask(:,:,:));
mask_P = squeeze(roiMask.RPMask(:,:,:));
mask_kid = mask_cor | mask_med | mask_CS | mask_P;
pixel_count1 = sum(mask_kid(:));
mask_kid1 = repmat((mask_kid),[1,1,1,151]);

mask_cor = squeeze(roiMask.LCMask(:,:,:));
%mask_cor = repmat((mask_cor),[1,1,1,151]);
mask_med = squeeze(roiMask.LMMask(:,:,:));
mask_CS = squeeze(roiMask.LCSMask(:,:,:));
mask_P = squeeze(roiMask.LPMask(:,:,:));
mask_kid2 = mask_cor | mask_med | mask_CS | mask_P;
pixel_count2 = sum(mask_kid2(:));
mask_kid2 = repmat((mask_kid2),[1,1,1,151]);
% mask_med = repmat((mask_med),[1,1,1,151]);
% mask_kid = mask_cor | mask_med ;
mask_aif = squeeze(roiMask.AMask(:,:,:));
mask_aif = repmat((mask_aif),[1,1,1,151]);
reconart1 = zeros(256,256,34,151);
reconart1(mask_aif >0) = (phantom60(mask_aif >0));
%x_aif = sum(reconart1(:))/pixel_count;
x_aif = squeeze(sum(sum(sum(abs(reconart1),1),2),3));
%x5_aif = smooth(x5_aif);x5_aif = smooth(x5_aif);x5_aif = smooth(x5_aif);
reconart1 = zeros(256,256,34,151);
reconart1(mask_kid1 >0) = (phantom60(mask_kid1 >0));
x_cor = squeeze(sum(sum(sum(abs(reconart1),1),2),3));


% load('mask_cortexPS3_kid1corr1.mat');
% mask_cortexkid1 = mask_cortex;
% pixel_count1 = sum(mask_cortexkid1(:));
% load('mask_cortexPS3_kid2corr1.mat');
% pixel_count2 = sum(mask_cortex(:));
% load('mask_aifPS3corr1.mat');
% load('/mnt/data/abhishek/PSdata/rdMeas/nav_recondcePS3.mat');

% parameters for 2 CM model
xFit(1) = 0.25;%45/(pixel_count1*1.32*1.32*6*1e-3);%0.5947;
xFit(2) = 0.44;
xFit(3) = 0.01;
XFit(4) = 0.01;
xFit(5) = 0.2;
xFit(6) = 4;
xFit(7) = 1;
%xFit(6) = 0.4;
% 

showPlot = 1;
range = 1:9;

TRes = 7*35*3.52e-3;%0.928;%mruDataset.TRes;
TR = 3.52e-3;%mruDataset.TR;
FA = 10*(pi/180);%mruDataset.FA;
clear phantom60;

% Get aorta and kidney signal
%increase the height

x_aif(20:end) = 1.8*x_aif(20:end);
x_aif(19) = (1.4)*x_aif(19);
x_aif(18) = 1.2*x_aif(18);
SI_aorta = x_aif;%averageFromMask(mruIm_cut,AMask)';
SI_kidney=x_cor;%averageFromMask(mruIm_cut,CMask)'; % Choose cortex for 3 compartment

% tkid = 2:3.8:217;
% tresnew = 1;
% tnew = 0.5:tresnew:216.5;
% tnew = reshape(tnew,length(tnew),1); %convert to column vector
% tkid = reshape(tkid,length(tkid),1);
%SI_kidney = interp1(tkid,SI_kidney,tnew,'pchip','extrap');

t = 0:(length(SI_aorta)-1);
t = t.*TRes;

T10_blood = 1.4;    % s
T10_kidney = 1.2;   % s
r1_blood = 4.5;     % (s^-1)*(mM^-1)
Htc_large = 0.41;

%// Find Cp_aorta
Cb = SI2C(SI_aorta,TR,FA,T10_blood,r1_blood);   % Concentration in blood
Cp_aorta = Cb/(1-Htc_large);% Concentration in plasma
save('Cp_aorta.mat','Cp_aorta');
%// Find Ckidney
Ckidney = SI2C(SI_kidney,TR,FA,T10_kidney,r1_blood);   % Concentration in blood
%save('Ckidney.mat','Ckidney');
Ckidney = Ckidney(1:110,:);
Cp_aorta = Cp_aorta(1:110,:);
t = t(1:110);
%[xFit,resnorm,xdata1,ydata] = FitThreeCompartment(Ckidney,Cp_aorta,t);
%[xFit,resnorm,xdata1,ydata] = Fit_modtofts(Ckidney,Cp_aorta,t);
figure

Cp_aorta = reshape(Cp_aorta,length(Cp_aorta),1);
xdata1 = [Cp_aorta,t'];
ydata = Ckidney;

% generate signals to be used
if (showPlot)
%     
    [Cfit, Cart, Ctub]= mod_sourbron(xFit,xdata1);
    plot(t,Cfit,'*k');
    hold on
    plot(t,Cart,'--r');
    plot(t,Ctub,'--b');
    %plot(t,ydata,'*k');
    plot(t,Cp_aorta,'+r');
    %plot(t,abs(Cfit-ydata),'+m');
    hold off
    legend('newParenchyma','IV','EV','AIF')
    title(['K21 = ' num2str(xFit(1)) ' ' 'vb = ' num2str(xFit(2)) ]);% ' resnorm = ' num2str(resnorm) ]);%' GFR =' num2str(round(pixel_count*xFit(1)*1.4*1.4*6*1e-3))])
    xlabel('Time (s)')
    ylabel('Concentration (mM)')
end

% converting cocentartion to SI
Cp_aorta = Cp_aorta*(1-Htc_large); 
SIaif = C2SI(Cp_aorta,TR,FA,T10_blood,r1_blood,SI_aorta(1,1));
save('Ckidney.mat','Cfit');
SIkidney1 = C2SI(Cfit,TR,FA,T10_kidney,r1_blood,SI_kidney(1,1));


xFit(1) = 0.25;%45/(pixel_count1*1.32*1.32*6*1e-3);%0.5947;
xFit(2) = 0.44;
xFit(3) = 0.01;
XFit(4) = 0.01;
xFit(5) = 0.2;
xFit(6) = 4;
xFit(7) = 1;
%xFit(6) = 0.4;
% 

showPlot = 1;
range = 1:9;

TRes = 7*35*3.52e-3;%0.928;%mruDataset.TRes;
TR = 3.52e-3;%mruDataset.TR;
FA = 10*(pi/180);%mruDataset.FA;
clear phantom60;

% Get aorta and kidney signal
SI_aorta = x_aif;%averageFromMask(mruIm_cut,AMask)';
SI_kidney=x_cor;%averageFromMask(mruIm_cut,CMask)'; % Choose cortex for 3 compartment

% tkid = 2:3.8:217;
% tresnew = 1;
% tnew = 0.5:tresnew:216.5;
% tnew = reshape(tnew,length(tnew),1); %convert to column vector
% tkid = reshape(tkid,length(tkid),1);
%SI_kidney = interp1(tkid,SI_kidney,tnew,'pchip','extrap');

t = 0:(length(SI_aorta)-1);
t = t.*TRes;

T10_blood = 1.4;    % s
T10_kidney = 1.2;   % s
r1_blood = 4.5;     % (s^-1)*(mM^-1)
Htc_large = 0.41;

%// Find Cp_aorta
Cb = SI2C(SI_aorta,TR,FA,T10_blood,r1_blood);   % Concentration in blood
Cp_aorta = Cb/(1-Htc_large);                    % Concentration in plasma

%// Find Ckidney
Ckidney = SI2C(SI_kidney,TR,FA,T10_kidney,r1_blood);   % Concentration in blood

Ckidney = Ckidney(1:110,:);
Cp_aorta = Cp_aorta(1:110,:);
t = t(1:110);
%[xFit,resnorm,xdata1,ydata] = FitThreeCompartment(Ckidney,Cp_aorta,t);
%[xFit,resnorm,xdata1,ydata] = Fit_modtofts(Ckidney,Cp_aorta,t);
figure

Cp_aorta = reshape(Cp_aorta,length(Cp_aorta),1);
xdata1 = [Cp_aorta,t'];
ydata = Ckidney;

% generate signals to be used
if (showPlot)
%     
    [Cfit, Cart, Ctub]= mod_sourbron(xFit,xdata1);
    plot(t,Cfit,'*k');
    hold on
    plot(t,Cart,'--r');
    plot(t,Ctub,'--b');
    %plot(t,ydata,'*k');
    plot(t,Cp_aorta,'+r');
    %plot(t,abs(Cfit-ydata),'+m');
    hold off
    legend('newParenchyma','IV','EV','AIF')
    title(['K21 = ' num2str(xFit(1)) ' ' 'vb = ' num2str(xFit(2)) ]);% ' resnorm = ' num2str(resnorm) ]);%' GFR =' num2str(round(pixel_count*xFit(1)*1.4*1.4*6*1e-3))])
    xlabel('Time (s)')
    ylabel('Concentration (mM)')
end

% converting cocentartion to SI
Cp_aorta = Cp_aorta*(1-Htc_large); 
%SIaif = C2SI(Cp_aorta,TR,FA,T10_blood,r1_blood,SI_aorta(1,1));
SIkidney2 = C2SI(Cfit,TR,FA,T10_kidney,r1_blood,SI_kidney(1,1));


load('mask_cortexPS3_kid1corr1.mat');
mask_cortexkid1 = mask_cortex;
pixel_count1 = sum(mask_cortexkid1(:));
load('mask_cortexPS3_kid2corr1.mat');
pixel_count2 = sum(mask_cortex(:));
load('mask_aifPS3corr1.mat');
load('/mnt/data/abhishek/PSdata/rdMeas/nav_recondcePS3.mat');

% put signal in there respective ROI's
% load('phantom60.mat');
% mask_aif = squeeze(roiMask.AMask(:,:,:));
% mask_cortex = squeeze(roiMask.RCMask(:,:,:));
% mask_medulla = squeeze(roiMask.RMMask(:,:,:));
% mask_CS = squeeze(roiMask.RCSMask(:,:,:));
% mask_P = squeeze(roiMask.RPMask(:,:,:));
% mask_kidney1 = mask_cortex | mask_medulla |mask_CS | mask_P ;
% mask_cortex = squeeze(roiMask.LCMask(:,:,:));
% mask_medulla = squeeze(roiMask.LMMask(:,:,:));
% mask_CS = squeeze(roiMask.LCSMask(:,:,:));
% mask_P = squeeze(roiMask.LPMask(:,:,:));
% mask_kidney2 = mask_cortex | mask_medulla |mask_CS | mask_P ;
% phantom_pre = squeeze(phantom60(:,:,:,1));
%phantom_pre = repmat(phantom_pre,[1,1,1,348]);
mask_aif = (mask_aif);
mask_kid1 = (mask_cortexkid1);
mask_kid2 = (mask_cortex);
%load('aorta_aif.mat');
%load('left_kid.mat');
SIaif_scale = ((SIaif./min(SIaif)));
SIkidney_scale1 = ((SIkidney1./min(SIkidney1)));
SIkidney_scale2 = ((SIkidney2./min(SIkidney2)));
%phantom_pre45 = double(repmat(insp_sos,[1,1,1,110]));
phantom_pre1 = TP1;%insp_sos(:,:,1:42);
phantom_pre = insp_sos(:,:,1:42);
% for loop to generate 110 time points
for i = 1:110
    
    phantom_pre1(squeeze(mask_aif(:,:,:))>0) = SIaif_scale(i).*phantom_pre(squeeze(mask_aif(:,:,:))>0);
    phantom_pre1(mask_kid1>0) = SIkidney_scale1(i).*phantom_pre(mask_kid1>0);
    phantom_pre1(mask_kid2>0) = SIkidney_scale2(i).*phantom_pre(mask_kid2>0);
    phantom_pre45(:,:,:,i) = phantom_pre1;
  
end;
% phantom_pre1 = phantom_pre;
% for i = 1:110
%     
%     %phantom_pre(mask_aif>0) = SIaif_scale(i).*phantom_pre(mask_aif>0);
%     phantom_pre1(mask_kidney>0) = SIkidney_scale(i).*phantom_pre(mask_kidney>0);
%     phantom_pre45(:,:,:,i) = phantom_pre1;
%   
% end;
% check if same signal is generated
% mask_cor = squeeze(roiMask.RCMask(:,:,:));
% %mask_cor = repmat((mask_cor),[1,1,1,151]);
% mask_med = squeeze(roiMask.RMMask(:,:,:));
% mask_CS = squeeze(roiMask.RCSMask(:,:,:));
% mask_P = squeeze(roiMask.RPMask(:,:,:));
% mask_kid = mask_cor | mask_med | mask_CS | mask_P;
% pixel_count = sum(mask_kid(:));
mask_kid = repmat((mask_kid1),[1,1,1,110]);
% mask_aif = squeeze(roiMask.AMask(:,:,:));
mask_aif = repmat((mask_aif),[1,1,1,110]);
reconart1 = zeros(288,288,42,110);
reconart1(mask_aif >0) = (phantom_pre45(mask_aif >0));
x5_aif = squeeze(sum(sum(sum(abs(reconart1),1),2),3));
%x5_aif = smooth(x5_aif);x5_aif = smooth(x5_aif);x5_aif = smooth(x5_aif);
reconart1 = zeros(288,288,42,110);
reconart1(mask_kid >0) = (phantom_pre45(mask_kid >0));
x5_cor = squeeze(sum(sum(sum(abs(reconart1),1),2),3));
% signal intensity curves
figure;
plot(x5_cor,'*k');
hold on
plot(x5_aif,'+r');

% check concentration curve
Cb = SI2C(x5_aif,TR,FA,T10_blood,r1_blood);   % Concentration in blood
Cp_aorta = Cb/(1-Htc_large);                    % Concentration in plasma

%// Find Ckidney
Ckidney = SI2C(x5_cor,TR,FA,T10_kidney,r1_blood);   % Concentration in blood

figure;
plot(Cp_aorta,'-g');
hold on;
plot(Ckidney,'-r');

save('/mnt/data/abhishek/phantom_pre45axial/phantom45_realistic2PS3.mat','phantom_pre45','mask_kid1','mask_aif','mask_kid2','-v7.3');
% clear all; 
% clc;