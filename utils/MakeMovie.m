function MakeMovie(imgt, filename,maxI)
if nargin<2
    filename = ['random_',randi(10^8)];
end
if nargin<3
    maxI=max(abs(imgt(:)))/1.05;
end
v = VideoWriter(filename);
open(v);
%Generate initial data and set axes and figure properties.

h=figure;

% axis tight manual 
% set(gca,'nextplot','replacechildren'); 
% waitfor(h)
%Create a set of frames and write each frame to the file.
nt = size(imgt,3);

for t = 1:nt
   imshow(abs(imgt(:,:,t)),[0,maxI]);
   frame = getframe;
   writeVideo(v,frame);
end
close(v);
close(h);
end