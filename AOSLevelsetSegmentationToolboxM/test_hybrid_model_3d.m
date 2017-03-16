%%
clear all; close all; clc
addpath('data'); 
load head_ct; %V; 

%%
% V = abs(imcnst);
V = impocs(:,:,1:48);
propagation_weight = 1e-3; 
GAC_weight = .02; 
% g = ones(size(V)); % linear diffusion 
g = ac_gradient_map(V,1,1); 
delta_t = 1; 
% mu = 1200; 
mu = 6.7e-5;

margin = 1; 
center = [129, 115, 24]; 
% center = [67,43, 25]; 
phi = zeros(size(V)); 
phi(center(1)-margin:center(1)+margin,...
    center(2)-margin:center(2)+margin,...
    center(3)-margin:center(3)+margin) = 1; 
%%
figure(100)
for i = 1:30
    phi = ac_hybrid_model(V-mu, phi-.5, propagation_weight, GAC_weight, g, ...
        delta_t, 1); 
    if exist('h','var') && all(ishandle(h)), delete(h); end
    iso = isosurface(phi,0);
    h = patch(iso,'edgecolor','r','facecolor','w');  axis equal;  view(3); 
    set(gcf,'name', sprintf('#iters = %d',i));
    drawnow; 
end

%%
figure;
% slice = [10,15,20,25,30,35,40,45];
 slice = (1:8)+32;
%  slice = 8:15;
for i = 1:8
    subplot(2,4,i); imshow(V(:,:,slice(i)),[]); hold on; 
%     contour(phi(:,:,slice(i)),[0,0],'r');
    c = contour(phi(:,:,slice(i)),[0,0]);
    [~,pt]=zy_plot_contours(c,'linewidth',2);
end