function [u] = prox_mask(g,sigma,mask)
%solve 0.5*||u-g||^2+mu*||C(R_b(u))||_* w.r.t u

%size(g) = nx,ny,T
%size(mask)=4; in cell format, for arota,left and right kidneys
[gx, gy, T]=size(g);
g = reshape(g,[],T);
N = length(mask);

for n=1:N
    inds = mask{n};
    [UU,SS,VV]=svd(g(inds,:),0);
    S  = diag(SS);
    mu = sigma*(sqrt(numel(inds))+sqrt(T));
    S  = (S-mu).*(S>mu);
    SS(1:length(S),1:length(S)) = diag(S);
    g(inds,:) = UU*SS*VV';
end
u=reshape(g,gx,gy,T);


