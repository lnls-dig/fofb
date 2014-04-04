function [num, den] = buildcic(decfactor, nstages, ndelays)

cic_coeff = tf(mfilt.cicdecim(decfactor, ndelays, nstages));
num = interp1(1:length(cic_coeff), cic_coeff, linspace(1, length(cic_coeff), nstages*ndelays+1));
num((ceil(length(num)/2)+1):end) = num((floor(length(num)/2)):-1:1);
num = num/sum(num);

den = [1 zeros(1,length(num)-1)];