function [sysd, fit, sys_array] = idcorrps(data, idmethod, order, idvalfactor, nk, Ts)

if nargin < 4 || isempty(idvalfactor)
    idvalfactor = 2/3;
end

if nargin < 5
    nk = [];
end

if nargin < 6
    Ts = -1;
end

ncorr = length(data);

sysd = zpk([],[],1,-1);
fit = zeros(ncorr,1);

for i=1:ncorr
    
    npts_ident = floor(size(data{i},1)*idvalfactor);

    data_ident = data{i}(1:npts_ident);
    data_val = data_ident; data{i}(npts_ident+1:end);
    
    if isempty(nk)
        nk_ = delayest(data_ident);
    else
        nk_ = nk;
    end
    sys = idmethod(data_ident, [order nk_]);
    sys_array{i} = sys;
    syszpk = zpk(sys);
    syszpk = syszpk(1,1);
    syszpk.Ts = Ts;
    syszpk = syszpk/dcgain(syszpk);
    sysd(i,i) = syszpk;
    dly = length(find(syszpk.p{1} == 0));
    syszpk.p = {syszpk.p{1}((syszpk.p{1} ~= 0))};
    [~, fit(i)] = compare(data_val, sys);
    fprintf('%s: %0.1f%% \n', data_ident.OutputName{1}, fit(i));
end