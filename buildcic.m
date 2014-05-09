function [num, den] = buildcic(decfactor, nstages, ndelays, method, density)

if nargin < 4 || isempty(method)
    method = 'freqz';
end
if nargin < 5 || isempty(density)
    density = 10;
end

cic_coeff = tf(mfilt.cicdecim(decfactor, ndelays, nstages));
cic_coeff = cic_coeff/sum(cic_coeff);
n = length(cic_coeff);
ndecim = nstages*ndelays;

if strcmpi(method, 'freqz')
    [h,w] = freqz(cic_coeff, [1 zeros(1,n-1)], density*ceil(n/decfactor)*decfactor);

    h = h(1:n/decfactor*density-1);
    w = linspace(0,pi,density*n/decfactor);
    w = w(1:end-1);

    num =invfreqz(h,w,ndecim,0);
    num = num/sum(num);
    den = [1 zeros(1,length(num)-1)];

elseif strcmpi(method, 'decim')
    num = interp1(1:length(cic_coeff), cic_coeff, linspace(1, length(cic_coeff), ndecim+1-rem(ndecim,2)));
    num((ceil((length(num)-1)/2)+1):end) = num((floor((length(num)-1)/2)+1):-1:1);
    num = num/sum(num);
    den = [1 zeros(1,length(num)-1)];

end