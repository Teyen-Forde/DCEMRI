function k=traj(read,views)
disp('design radial trajectory')
k = complex(zeros(read,views),zeros(read,views));
% traj=linspace(-1/2,1/2,read);
traj=(-read/2+0.5:1:read/2-0.5)/read;
deg_skip = (sqrt(5)-1)/2*pi;
theta = (0:views-1)*deg_skip;

for i=1:views
    for j=1:read
        kx = sin( theta(i) ) * traj(j);
        ky = cos( theta(i) ) * traj(j);
        k(j,i) = complex(kx,ky);
    end
end
disp('done.')
end