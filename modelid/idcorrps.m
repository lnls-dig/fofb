function [sysd, fit] = idcorrps(indata, outdata, names, idmethod, order, percent, nk, Ts)

if nargin < 6 || isempty(percent)
    percent = 80;
end

if nargin < 7
    nk = [];
end

if nargin < 8
    Ts = -1;
end

ncorr = size(indata,2);
npts_ident = floor(size(indata,1)*percent/100);

idin = indata(1:npts_ident,:);
idout = outdata(1:npts_ident,:);
valin = indata(npts_ident+1:end,:);
valout = outdata(npts_ident+1:end,:);

sysd = zpk([],[],1,-1);
fit = zeros(ncorr,1);

for i=1:size(indata,2)
    if isempty(nk)
        nk_ = delayest(iddata(outdata(:,i), indata(:,i)));
    else
        nk_ = nk;
    end
    sys = idmethod(iddata(idout(:,i), idin(:,i)), [order nk_]);
    syszpk = zpk(sys);
    syszpk = syszpk(1,1);
    syszpk.Ts = Ts;
    syszpk = syszpk/dcgain(syszpk);
    sysd(i,i) = syszpk;
    dly = length(find(syszpk.p{1} == 0));
    syszpk.p = {syszpk.p{1}((syszpk.p{1} ~= 0))};
    syszpk = d2c(zpk(syszpk));
    syszpk.inputdelay = Ts*dly;
    [~, fit(i)] = compare(iddata(valout(:,i), valin(:,i)), sys);
    fprintf('%s: %0.1f%% \n', names{i}, fit(i));
end