addpath('utils')
load kidney_phantomgeneration/phantom60
slice = squeeze(phantom60(:,:,17,1:2:end));
clear phantom60
[noise,sigma]=awgnc(slice,20);
data = single(slice+noise);
%% LLR -sliding
T=size(data,3);
b=16;
mu = sigma*(b+sqrt(T));
tic
[slice_LLRs,sv_LLRs] = prox_LLR(data,mu,[b,b]);%time 240s
toc
%% GLR
T=size(data,3);
b=size(data,1);
mu = sigma*(b+sqrt(T));
tic
[slice_GLR,sv_GLR] = prox_LLR(data,mu,[b,b]);%time 2.5s
toc
%% LLR -distinct -randshift
T=size(data,3);
b=16;
mu = sigma*(b+sqrt(T));
tic
[slice_LLRd,svals] = llr_thresh(data,mu,[b,b]);%time 1.3s
sv_LLRd=squeeze(svals(3,6,:));
toc
%%
r=1;
for f=1:76
    tmp=isame(slice(:,:,f),slice_LLRd(:,:,f));
    if tmp<r
        r=tmp;
        frame=f;
    end
end
%%
% frame=51;
data_frame = data(:,:,frame);
LLR_frame = slice_LLRd(:,:,frame);
GLR_frame = slice_GLR(:,:,frame);
figure,imshow(abs([slice(:,:,frame),data_frame,GLR_frame,LLR_frame;...
                   20*(slice(:,:,frame)-slice(:,:,frame)),...
                   20*(data_frame-slice(:,:,frame)),...
                   20*(GLR_frame-slice(:,:,frame)),...
                   20*(LLR_frame-slice(:,:,frame))]),[])
ftitle('Ground truth, noisy image, GLR recon, LLR recon')
xlabel('20x difference'); faxis;
r1=isame(slice(:,:,frame),data_frame);
r2=isame(slice(:,:,frame),GLR_frame);
r3=isame(slice(:,:,frame),LLR_frame);

%% output figure
figure,plot(log10(sv_GLR./max(sv_GLR)),'LineWidth',3, 'Color','k')
hold on
plot(log10(sv_LLRd./max(sv_LLRd)),'LineWidth',3,'Color','g')
legend('GLR','LLR');faxis;
ftitle('normalized singular values in log scale')