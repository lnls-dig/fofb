function [a, b, fit] = idsiso(datae, datav, order, sysarray, fixedsys)

na = order(1);
nb = order(2);
nk = order(3);

nu = length(datae);
a = zeros(nu, max(na+1,nb+nk));
b = zeros(nu, max(na+1,nb+nk));
fit = zeros(nu,1);
for i = 1:nu
    if iscell(sysarray)
        sys = sysarray{i};
    else
        sys = sysarray;
    end
    [datae_, datav_] = idpreprocessdata(detrend(datae{i},0), detrend(datav{i},0), 1, 1, sys, fixedsys, size(datae{i},1));
    m = arx(datae_, order, 'Focus', 'Simulation');
    [~, fit(i)] = compare(m, datav_, Inf);
    mtf = tf(m, 'measured');
    den = mtf.den{1};
    num = mtf.num{1};
    a(i,:) = den;
    b(i,:) = num;
end