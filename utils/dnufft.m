function FT = nufft(k,imSize,USE,SCALE)
if nargin < 4
    SCALE=false;
end
N=[imSize,imSize];
dims = size(k);
disp('design nufft operator')
if USE == 'c'%use cpu
    if SCALE
        w = prod(N)/prod(dims)*pi/2;
    else
        w = 1;
    end
    FT = NUFFT(k,w,1,0,N,2);
elseif USE == 'g'
    if SCALE
        w = ones(dims)*(prod(N)/prod(dims)*pi);
    else
        w = ones(dims);
    end
    FT = NUFFT_GPU(N,N*2, 5.5, dims(2),k,w);
else
    error('only have cpu or gpu mode!')
end
disp('done.')
end