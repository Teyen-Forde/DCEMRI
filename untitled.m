load imgt27;
load outCSVt2
load outCSVt4
load outCSVt6
load outCSV4
frame=21;
%% figure 1 and 2 conclusion: training is OK, at least as good as using exact basis; truncated at 4
figure(1)%obj error convergence of Vt2,4,6 and V4
plot(outCSVt2.time,outCSVt2.objerr,'LineWidth',2); hold on
plot(outCSVt4.time,outCSVt4.objerr,'LineWidth',2); hold on
plot(outCSVt6.time,outCSVt6.objerr,'LineWidth',2); hold on
plot(outCSV4.time,outCSV4.objerr,'LineWidth',2); hold off
title('obj err |Ax-b|')
xlabel('t/s')
legend('TB2','TB4','TB6','EB4')
figure(2)%recons of Vt2,4,6,and V4
%frame frame is the maximum
BIG=[];
BIG=[BIG,abs([outCSVt2.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSVt2.CS(:,:,frame))])];
BIG=[BIG,abs([outCSVt4.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSVt4.CS(:,:,frame))])];
BIG1=BIG;
BIG=[];
BIG=[BIG,abs([outCSVt6.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSVt6.CS(:,:,frame))])];
BIG=[BIG,abs([outCSV4.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSV4.CS(:,:,frame))])];
BIG2=BIG;
BIG=[BIG1;BIG2];
figure,imshow(imgt(:,:,frame),[]); title('Frame 21')
ma=max(max(imgt(:,:,frame)));
figure(2), imshow(BIG,[0,ma]);
%%
load outCS 
load outCSV4 
load outCSVt4 
load outLSI 
load outLSA
figure(3)% comparing CS, ktpca, lsi, lsa
plot(outCS.time,outCSVt2.objerr,'LineWidth',2); hold on
plot(outCSVt4.time,outCSVt4.objerr,'LineWidth',2); hold on
plot(outLSI.time,outLSI.objerr,'LineWidth',2); hold on
plot(outLSA.time,outLSA.objerr,'LineWidth',2); hold off
title('obj err |Ax-b|')
xlabel('t/s')
legend('CS','kt-PCA','LSI','LSA')

%%
%frame frame is the maximum
BIG=[];
BIG=[BIG,abs([outCS.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCS.CS(:,:,frame))])];
BIG=[BIG,abs([outCSVt4.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSVt4.CS(:,:,frame))])];
BIG1=BIG;
BIG=[];
tmp1=outLSI.L(:,:,frame)+outLSI.S(:,:,frame);
BIG=[BIG,abs([tmp1,5*abs(imgt(:,:,frame)-tmp1)])];
tmp2=outLSA.L(:,:,frame)+outLSA.S(:,:,frame);
BIG=[BIG,abs([tmp2,5*abs(imgt(:,:,frame)-tmp2)])];
BIG2=BIG;
BIG=[BIG1;BIG2];
ma=max(max(imgt(:,:,frame)));
figure(4), imshow(BIG,[0,ma]);colormap jet
figure,imshow(imgt(:,:,frame),[]); title('Frame frame'), colormap jet

%% Signal intensity 
addpath('kidney_phantomgeneration/')
load('mask_aifPS3corr1.mat')
mask_aif=mask_aif(:,:,27);
load('mask_cortexPS3_kid1corr1.mat')
mask_kid1=mask_cortex(:,:,27);
load('mask_cortexPS3_kid2corr1.mat')
mask_kid2=mask_cortex(:,:,27);

load EhnPht
EhnPht = squeeze(EhnPht(:,:,27,:));
%%
inds_aif=find(mask_aif);
inds_kid1=find(mask_kid1);
inds_kid2=find(mask_kid2);
phant=abs(reshape(EhnPht,[],nt));
SI_aif1=sum(phant(inds_aif,:));
SI_1kid1=sum(phant(inds_kid1,:));
SI_2kid1=sum(phant(inds_kid2,:));

tmp2 = outCS.CS;
% tmp2 = CS;
phant=abs(reshape(rot90(tmp2,-1),[],nt));
SI_aif2=sum(phant(inds_aif,:));
SI_1kid2=sum(phant(inds_kid1,:));
SI_2kid2=sum(phant(inds_kid2,:));

tmp3 = outCSV4.CS;
phant=abs(reshape(rot90(tmp3,-1),[],nt));
SI_aif3=sum(phant(inds_aif,:));
SI_1kid3=sum(phant(inds_kid1,:));
SI_2kid3=sum(phant(inds_kid2,:));

tmp4= outCSVt4.CS;
% tmp4 = CS;
phant=abs(reshape(rot90(tmp4,-1),[],nt));
SI_aif4=sum(phant(inds_aif,:));
SI_1kid4=sum(phant(inds_kid1,:));
SI_2kid4=sum(phant(inds_kid2,:));

tmp5=outLSA.L+outLSA.S;
phant=abs(reshape(rot90(tmp5,-1),[],nt));
SI_aif5=sum(phant(inds_aif,:));
SI_1kid5=sum(phant(inds_kid1,:));
SI_2kid5=sum(phant(inds_kid2,:));

% tmp6=outLSI.L+outLSI.S;
tmp6 = LplusS;
phant=abs(reshape(rot90(tmp6,-1),[],nt));
SI_aif6=sum(phant(inds_aif,:));
SI_1kid6=sum(phant(inds_kid1,:));
SI_2kid6=sum(phant(inds_kid2,:));
%%

figure(5);
plot(SI_aif1,'LineWidth',2); hold on; 
plot(SI_aif2,'LineWidth',2);plot(SI_aif3,'LineWidth',2);
plot(SI_aif4,'LineWidth',2);plot(SI_aif5,'LineWidth',2); 
plot(SI_aif6,'LineWidth',2); hold off
title('aif signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA','LSI')

figure(6),
subplot(221)
plot(SI_1kid1,'LineWidth',2), hold on;
plot(SI_1kid2,'LineWidth',2);plot(SI_1kid3,'LineWidth',2);
plot(SI_1kid4,'LineWidth',2);plot(SI_1kid5,'LineWidth',2);
plot(SI_1kid6,'LineWidth',2);hold off;
title('kid1 signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA','LSI')

subplot(222)
plot(SI_1kid1,'LineWidth',2), hold on;
plot(SI_1kid2,'LineWidth',2);plot(SI_1kid3,'LineWidth',2);
plot(SI_1kid4,'LineWidth',2);plot(SI_1kid5,'LineWidth',2);
plot(SI_1kid6,'LineWidth',2);hold off;

subplot(223)
plot(SI_2kid1,'LineWidth',2), hold on;
plot(SI_2kid2,'LineWidth',2);plot(SI_2kid3,'LineWidth',2);
plot(SI_2kid4,'LineWidth',2);plot(SI_2kid5,'LineWidth',2);
plot(SI_2kid6,'LineWidth',2);hold off;
title('kid2 signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA','LSI')

subplot(224)
plot(SI_2kid1,'LineWidth',2), hold on;
plot(SI_2kid2,'LineWidth',2);plot(SI_2kid3,'LineWidth',2);
plot(SI_2kid4,'LineWidth',2);plot(SI_2kid5,'LineWidth',2);
plot(SI_2kid6,'LineWidth',2);hold off;

%% d
t=1:18:76;
LplusSd=[LplusS(:,:,t(1)),LplusS(:,:,t(2)),LplusS(:,:,t(3)),LplusS(:,:,t(4)),LplusS(:,:,t(5))];
Ld=[L(:,:,t(1)),L(:,:,t(2)),L(:,:,t(3)),L(:,:,t(4)),L(:,:,t(5))];
Sd=[S(:,:,t(1)),S(:,:,t(2)),S(:,:,t(3)),S(:,:,t(4)),S(:,:,t(5))];
TSd=[S(:,:,t(1)+1)-S(:,:,t(1)),S(:,:,t(2)+1)-S(:,:,t(2)),S(:,:,t(3)+1)-S(:,:,t(3)),S(:,:,t(4)+1)-S(:,:,t(4)),S(:,:,t(5)+1)-S(:,:,t(5))];


figure(8);
subplot(4,1,1),imshow(abs(LplusSd),[]);ylabel('L+S')
subplot(4,1,2),imshow(abs(Ld),[]);ylabel('L')
subplot(4,1,3),imshow(abs(Sd),[]);ylabel('S')
subplot(4,1,4),imshow(abs(TSd),[]);ylabel('TS')

%% noise level 30dB
load outCS30 
load outCSV30 
load outCSVt30 
load outLSA30

load outCS20 
load outCSV20 
load outCSVt20 
load outLSA20
inds_aif=find(mask_aif);
inds_kid1=find(mask_kid1);
inds_kid2=find(mask_kid2);
phant=abs(reshape(EhnPht,[],nt));
SI_aif1=sum(phant(inds_aif,:));
SI_1kid1=sum(phant(inds_kid1,:));
SI_2kid1=sum(phant(inds_kid2,:));

tmp2 = outCS20.CS;
phant=abs(reshape(rot90(tmp2,-1),[],nt));
SI_aif2=sum(phant(inds_aif,:));
SI_1kid2=sum(phant(inds_kid1,:));
SI_2kid2=sum(phant(inds_kid2,:));

tmp3 = outCSV20.CS;
phant=abs(reshape(rot90(tmp3,-1),[],nt));
SI_aif3=sum(phant(inds_aif,:));
SI_1kid3=sum(phant(inds_kid1,:));
SI_2kid3=sum(phant(inds_kid2,:));

tmp4= outCSVt20.CS;
phant=abs(reshape(rot90(tmp4,-1),[],nt));
SI_aif4=sum(phant(inds_aif,:));
SI_1kid4=sum(phant(inds_kid1,:));
SI_2kid4=sum(phant(inds_kid2,:));

tmp5=outLSA20.L+outLSA20.S;
phant=abs(reshape(rot90(tmp5,-1),[],nt));
SI_aif5=sum(phant(inds_aif,:));
SI_1kid5=sum(phant(inds_kid1,:));
SI_2kid5=sum(phant(inds_kid2,:));


figure(5);
plot(SI_aif1,'LineWidth',2); hold on; 
plot(SI_aif2,'LineWidth',2);plot(SI_aif3,'LineWidth',2);
plot(SI_aif4,'LineWidth',2);plot(SI_aif5,'LineWidth',2); hold off
title('aif signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA')

figure(6),
subplot(221)
plot(SI_1kid1,'LineWidth',2), hold on;
plot(SI_1kid2,'LineWidth',2);plot(SI_1kid3,'LineWidth',2);
plot(SI_1kid4,'LineWidth',2);plot(SI_1kid5,'LineWidth',2);hold off;
title('kid1 signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA')

subplot(222)
plot(SI_1kid1,'LineWidth',2), hold on;
plot(SI_1kid2,'LineWidth',2);plot(SI_1kid3,'LineWidth',2);
plot(SI_1kid4,'LineWidth',2);plot(SI_1kid5,'LineWidth',2);hold off;

subplot(223)
plot(SI_2kid1,'LineWidth',2), hold on;
plot(SI_2kid2,'LineWidth',2);plot(SI_2kid3,'LineWidth',2);
plot(SI_2kid4,'LineWidth',2);plot(SI_2kid5,'LineWidth',2);hold off;
title('kid2 signal intensity')
legend('GroundTruth','CS', 'kt-PCA\_exact','kt-PCA','LSA')

subplot(224)
plot(SI_2kid1,'LineWidth',2), hold on;
plot(SI_2kid2,'LineWidth',2);plot(SI_2kid3,'LineWidth',2);
plot(SI_2kid4,'LineWidth',2);plot(SI_2kid5,'LineWidth',2);hold off;



%% noise level 20dB

relerr=zeros(1,nt);
for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),tmp2(:,:,tt));
end

figure(9),plot(relerr,'LineWidth',2); hold on

for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),tmp3(:,:,tt));
end
plot(relerr,'LineWidth',2);

for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),tmp4(:,:,tt));
end
plot(relerr,'LineWidth',2);

for tt=1:nt
    relerr(tt)=isame(imgt(:,:,tt),tmp5(:,:,tt));
end
plot(relerr,'LineWidth',2); hold off;
legend('CS', 'kt-PCA\_exact','kt-PCA','LSA')
title('Frame rel err @ 20dB')




%%

BIG=[];
BIG=[BIG,abs([outCS30.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCS30.CS(:,:,frame))])];
BIG=[BIG,abs([outCSV30.CS(:,:,frame),5*abs(imgt(:,:,frame)-outCSV30.CS(:,:,frame))])];
BIG1=BIG;
BIG=[];
tmp1=outCSVt30.CS(:,:,frame);
BIG=[BIG,abs([tmp1,5*abs(imgt(:,:,frame)-tmp1)])];
tmp2=outLSA30.L(:,:,frame)+outLSA30.S(:,:,frame);
BIG=[BIG,abs([tmp2,5*abs(imgt(:,:,frame)-tmp2)])];
BIG2=BIG;
BIG=[BIG1;BIG2];
ma=max(max(imgt(:,:,frame)));
figure(4), imshow(BIG,[0,ma]);colormap jet
figure,imshow(imgt(:,:,frame),[]); title('Frame frame'), colormap jet
ylabel('30dB')













