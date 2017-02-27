function [u,svals] = prox_LLR(g,mu,block_dim)
%solve 0.5*||u-g||^2+mu*||C(R_b(u))||_* w.r.t u
%C : transform into Casorati form B^2 x T
%Rb : extract the B*B*T block
%mu : estimated by sigma*(sqrt(B)+sqrt(T)), sigma is noise standard
%deviation of image serials

%size(g) = nx*ny*T

if nargin < 3
    block_dim=[10,10];
end
[gx, gy, T]=size(g);
% mx=gx+bx-1;my=gy+by-1;
% gpad = zpad(g,[mx,my,T]);
A=im2row(g,block_dim);
N=size(A,1);%# of blocks
A=permute(A,[2,3,1]);
r=randi(N,1);
% s_vals = zeros(K, L);
for n=1:N
    [UU, SS, VV] = svd(A(:,:,n),0);
    S=diag(SS);
    if n==r
        svals=S;
    end
    S=(S-mu).*(S>mu);
    A(:,:,n) = UU*diag(S)*VV';
    
end
A=permute(A,[3,1,2]);
u = row2im(A,[gx,gy],block_dim);

