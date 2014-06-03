function [sysset, fit, noisesysset] = idcorrps(indata, outdata, names, idmethod, order, percent, nk)

if nargin < 6 || isempty(percent)
    percent = 80;
end

if nargin < 7
    nk = [];
end

ncorr = size(indata,2);
npts_ident = floor(size(indata,1)*percent/100);

idin = indata(1:npts_ident,:);
idout = outdata(1:npts_ident,:);
valin = indata(npts_ident+1:end,:);
valout = outdata(npts_ident+1:end,:);

sysset = ss([],[],[],[]);
noisesysset = ss([],[],[],[]);
fit = zeros(ncorr,1);

for i=1:size(indata,2)
    if isempty(nk)
        nk_ = delayest(iddata(outdata(:,i), indata(:,i)));
    end
    sys = idmethod(iddata(idout(:,i), idin(:,i)), [order nk_]);
    syszpk = zpk(sys);
    sysset(i,i) = syszpk(1,1);
    noisesysset(i,i) = syszpk(1,2);
    [~, fit(i)] = compare(iddata(valout(:,i), valin(:,i)), sys);
    fprintf('%s: %0.1f%% \n', names{i}, fit(i));
end