function res = mtimes(a,b)

if a.adjoint
    res = adjDx(b(:,:,:,1))+adjDy(b(:,:,:,2))+adjDt(b(:,:,:,3));
else
    % input 2d dynamic image
    [nx,ny,nt] = size(b);
    res = zeros(nx,ny,nt,3);
    res(:,:,:,1) = Dx(b);
    res(:,:,:,2) = Dy(b);
    res(:,:,:,3) = Dt(b);
end

function res = Dx(x)
res = x([2:end,end],:,:) - x;

function res = Dy(x)
res = x(:,[2:end,end],:) - x;

function res = Dt(x)
res = x(:,:,[2:end,end]) - x;


function res = adjDx(x)
res= x([1,1:end-1],:,:) - x;
res(1,:,:) = -x(1,:,:);
res(end,:,:) = x(end-1,:,:);

function res = adjDy(x)
res= x(:,[1,1:end-1],:) - x;
res(:,1,:) = -x(:,1,:);
res(:,end,:) = x(:,end-1,:);

function res = adjDt(x)
res= x(:,:,[1,1:end-1]) - x;
res(:,:,1) = -x(:,:,1);
res(:,:,end) = x(:,:,end-1);