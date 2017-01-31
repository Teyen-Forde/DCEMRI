function imshow_obo(imgs)
%imshow one by one
coils=size(imgs,3);
imgs=abs(imgs);

for nc=1: coils
    figure(100), imshow(imgs(:,:,nc),[]);
    title(num2str(nc));
      pause;
end
end