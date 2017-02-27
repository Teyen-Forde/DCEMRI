function res = meanfilter(Sig,len)
%len is the length of filter
if nargin < 2
    len = 5;
end
if mod(len,2)==0
    len = len+1;
end
ls = length(Sig);
if ls < len
    error('Signal length is too short to do mean filtering')
end

%pading signal
S(1:(len-1)/2)=Sig((len-1)/2:-1:1);
S((len-1)/2+1:(len-1)/2+ls)=Sig;
S((len-1)/2+ls+1:len+ls)=Sig(end:-1:end-(len-1)/2);
res = zeros(size(Sig));
for l=1:ls
    res(l) = mean(S(l:l+len));
end
% S = zeros(ls+len,1);
% h = ones(len,1)/len;
% res = conv(Sig,h,'valid');
