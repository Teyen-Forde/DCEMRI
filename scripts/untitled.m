nt=100;
FT=cell(nt,1);
imSize=256;
N=[imSize,imSize];
nspokes=21;
coils=24;
for tt=1:nt
    
    kt=k(:,(tt-1)*nspokes+1:tt*nspokes);
%     wt=abs(kt)./max(abs(kt(:)));
%         FT{tt}=NUFFT(kt,1,1,0,N,2);
    FT{tt}=NUFFT_GPU(N,N*2, 5.5, nspokes,kt,ones(size(kt)));
end


lowres_img = bart('nufft -i -d24:24:1 -t', traj, ksp);


lowres_ksp = bart('fft -u 7', lowres_img);

% zeropad to full size
ksp_zerop = bart('resize -c 0 256 1 256', lowres_ksp);

% ESPIRiT calibration
sens = bart('ecalib -m1', ksp_zerop);

img = bart('nufft -i  -t', traj, ksp);
ksp = bart('fft -u 7',img);
sens1 = bart('ecalib -k 4 -r 24 -m1', ksp);


deg_skip = (sqrt(5)-1)/2*pi;
k = zeros(read,views);
traj = linspace(-1/2, 1/2,read);
theta = (0:views-1)*deg_skip;
for i=1:views
    for j=1:read
        k(j,i) = complex(traj(j)*sin(theta(i)),traj(j)*cos(theta(i)));
    end
end

w = abs(k);

FF=NUFFT(k,1,1,0,N,2);

imcoils= FF'*(ksp.*repmat(w,[1,1,6]));

figure,imshow(sos(imcoils))

tic
[PCCp,outCS] = CS_tv(params,kdata,sens);
toc

figure
for i=1:100
    imshow(rot90(abs(PCCp(:,:,i)),2),[0 2.5e-3]);
    pause(0.1)
end