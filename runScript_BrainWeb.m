addpath('./utils')
addpath('./data')
% example
[x,y] = meshgrid([-64:63]/128);
k = x + i*y;
w = 1;
phase = 1;
imSize = [128,128];
shift = [0,0];
FT = NUFFT(k,w,phase,shift,imSize,2);
clear
lambda=[5e-3, 0.01,0.01,0.015];
% reqSNR=[0,35,30,25];
load BrainWebData
sigma=0.0145;
mu=[0.5,1,2,4];
%%
for count=1:4
%     if(reqSNR(count))
        [~,sigma]=awgnc(kdata(:,:,:,:,1),48);
%     else
%         noise = 0;
%         sigma = 1e-4;
%     end
%     noise_std=sqrt(noisePower/2); 
    rand('seed',0);
    rnoise=randn(size(kdata));
    rand('seed',1);
    inoise=randn(size(kdata));
    noise=sigma*complex(rnoise,inoise);
    data=kdata+noise;
    %% nlcg
    clear params
    load params_NLCG
    params.smaps = smaps;
    params.data = kdata;
    params.cgiter =300;
    params.TVspatWeight = 5e-4;
    params.dbout = 1;
    tic;
    [PCCn, out_NLCG,time_NLCG,numA] = simple_repcom_nlcg(params,basis);
    [RMSE_3PCn, RMSE_T2n,ntt2, TE_IMGn]=postprocess(out_NLCG,gtt2,gtt3,basis,echotimes,mask);
    toc;
    %% padm
    clear params
    load params_PADM
    params.smaps = squeeze(smaps);
    params.data = squeeze(data);
    params.basis = basis;
    params.maxiter =30;
    params.epsilon = sqrt(numel(kdata)+8*sqrt(numel(kdata)))*(sigma);
    params.lbd = 1.2e-5/1.5;

    x0 = zeros(numel(mask),3);
    tic
    [PCCp, out_PADM,time_PADM] = padm_CG(x0,params);
    toc
    %post processing
    [RMSE_3PCp, RMSE_T2p,ptt2,TE_IMGp]=postprocess(out_PADM,gtt2,gtt3,basis,echotimes,mask);
    
    %% plotting
    figure(1),
    semilogy(time_NLCG,RMSE_3PCn,'r','Linewidth',1.5); hold on
    semilogy(time_PADM,RMSE_3PCp,'b','Linewidth',1.5);
    legend('NLCG','ADMM');
    semilogy(time_NLCG(50:50:end),RMSE_3PCn(50:50:end),'ro','Linewidth',1.5)
    semilogy(time_PADM(50:50:end),RMSE_3PCp(50:50:end),'bo','Linewidth',1.5); hold off
    title('RelErr of 3PCs in real time','FontSize',12,'FontWeight','bold')
        xlabel('time/s','FontSize',12,'FontWeight','bold')

    figure(2),
    semilogy(time_NLCG,RMSE_T2n,'r','Linewidth',1.5); hold on
    semilogy(time_PADM,RMSE_T2p,'b','Linewidth',1.5);
    legend('NLCG','ADMM');
    semilogy(time_NLCG(50:50:end),RMSE_T2n(50:50:end),'ro','Linewidth',1.5)
    semilogy(time_PADM(50:50:end),RMSE_T2p(50:50:end),'bo','Linewidth',1.5); hold off
    title('RelErr of T2 in real time','FontSize',12,'FontWeight','bold')
    xlabel('time/s','FontSize',12,'FontWeight','bold')
    figure(3),
    subplot(211),imshow(abs([gtt2,ntt2,ptt2]),[]); colormap jet; set(gca,'CLim',[0,300]); colorbar;
    subplot(212),imshow(abs([gtt2-gtt2,ntt2-gtt2,ptt2-gtt2]),[]); colormap jet; set(gca,'CLim',[0,100]); colorbar;
    %     title('T2 and error maps: Groud Truth, NLCG, PADM')
    figure(4),imshow(abs([TE_IMGn(:,:,16),TE_IMGp(:,:,16)])),colormap jet
    title('TE16 of NLCG and PADM recons')
    figure(5),imshow3(abs(reshape(TE_IMGp,imsteps,imsteps,ETL)))
    
    %% saving
    save(['SNR48',num2str(mu(count))],'time_NLCG','RMSE_3PCn','RMSE_T2n','ntt2','TE_IMGn','time_PADM','RMSE_3PCp','ptt2','RMSE_T2p','TE_IMGp');
    %     save(['SNRPADM1',num2str(reqSNR(count))],'time_PADM','RMSE_3PCp','ptt2','ptt2');
end

%%
figure(1),
    subplot(211),imshow(abs([gtt2,ntt2,ptt2]),[]); colormap jet; set(gca,'CLim',[0,300]); colorbar;
    subplot(212),imshow(abs([gtt2-gtt2,ntt2-gtt2,ptt2-gtt2]),[]); colormap jet; set(gca,'CLim',[0,100]); colorbar;
