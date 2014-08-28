function [a, b, fit, rcorr] = idsiso(datae, datav, idmethod, order, sysarray, fixedsys)

nlags = 25;
nu = length(datae);

% na = order(1);
% nb = order(2);
% nk = order(3);
% a = zeros(nu, max(na+1,nb+nk));
% b = zeros(nu, max(na+1,nb+nk));

a = {};
b = {};

fit = zeros(nu,1);
rcorr = zeros(2*nlags+1, nu);
for i = 1:nu
    if iscell(sysarray)
        sys = sysarray{i};
    else
        sys = sysarray;
    end
    [datae_, datav_] = idpreprocessdata(detrend(datae{i},0), detrend(datav{i},0), 1, 1, sys, fixedsys);
    
    if isequal(idmethod, @impulse)
        m = idmethod(datae_, 'PW', order(1));
    elseif  isequal(idmethod, @n4sid)
        %nk = delayest(datae_);
        %nk = max(nk, 1);
        m = idmethod(datae_, 'best');%, 'nk', nk);
        %m = pem(m,datae_);
        %plot(fft(datae_)); pause
    else
        m = idmethod(datae_, order, 'Focus', 'Simulation');
    end
        
    [~, fit(i)] = compare(m, datav_, Inf);
    
    r = resid(m, datav_);
    rcorr(:,i) = xcorr(r.OutputData, r.InputData, nlags, 'coeff');
    
    mtf = tf(m, 'measured');
    den = mtf.den{1};
    num = mtf.num{1};
    
    a{i} = den;
    b{i} = num;
end