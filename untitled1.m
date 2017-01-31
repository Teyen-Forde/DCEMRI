for tt =1: nt
    img(:,:,:,tt) = FT{tt}'*(kdata(end/2-2:end/2+1,:,:,tt).*mask);
end
M=reshape(img,[],nt);
[Ut,St,Vt]=svd(M,0);
figure,plot(diag(St))

%%
M=reshape(EhnPht,[],nt);
[Ut,St,Vt]=svd(M,0);