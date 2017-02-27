function mask = drawROI(im,verbose)
if nargin < 2
    verbose=0;
end
figure(100), imshow(im,[]);
set(gcf, 'Position', get(0, 'Screensize'));
title('Select ROI');
h = imfreehand;
wait(h);
mask = createMask(h);
if verbose
    hold on;
    imshow(im.*~mask,[])
    hold off;
end
end