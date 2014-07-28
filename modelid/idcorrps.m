function [sysd, fit, sys_array] = idcorrps(data_id, data_val, idmethod, order, nk, Ts)

if nargin < 5
    nk = [];
end

if nargin < 6
    Ts = -1;
end

ncorr = length(data_id);

sysd = zpk([],[],1,-1);
fit = zeros(ncorr,1);

for i=1:ncorr
    if isempty(nk)
        nk_ = delayest(data_id{i});
    else
        nk_ = nk;
    end
    sys = idmethod(detrend(data_id{i},0), [order nk_]);
    sys_array{i} = sys;
    syszpk = zpk(sys);
    syszpk = syszpk(1,1);
    syszpk.Ts = Ts;
    syszpk = syszpk/dcgain(syszpk);
    sysd(i,i) = syszpk;
    dly = length(find(syszpk.p{1} == 0));
    syszpk.p = {syszpk.p{1}((syszpk.p{1} ~= 0))};
    [~, fit(i)] = compare(detrend(data_val{i},0), sys, Inf);
    fprintf('%s: %0.1f%% \n', data_id{i}.OutputName{1}, fit(i));
end