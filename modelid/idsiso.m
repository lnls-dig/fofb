function [models, fit, rcorr] = idsiso(datae, datav, idmethod, order, input_dynamics, commmon_dynamics)

nlags = 25;
nu = length(datae);

a = {};
b = {};

fit = zeros(nu,1);
rcorr = zeros(2*nlags+1, nu);
for i = 1:nu
    if iscell(input_dynamics)
        sys = input_dynamics{i};
    else
        sys = input_dynamics;
    end
    [datae_, datav_] = idpreprocessdata(detrend(datae{i},0), detrend(datav{i},0), 1, 1, sys, commmon_dynamics);
    
    if isequal(idmethod, @impulse)
        m = idmethod(datae_, 'PW', order(1));
    elseif  isequal(idmethod, @n4sid)
        nk = delayest(datae_);
        nk = max(nk, 1);
        m = idmethod(datae_, 'best', 'nk', nk);
    else
        m = idmethod(datae_, order, 'Focus', 'Simulation');
    end

    % Fit quality
    [~, fit(i)] = compare(m, datav_, Inf);

    % Residual analysis
    r = resid(m, datav_);
    rcorr(:,i) = xcorr(r.OutputData, r.InputData, nlags, 'coeff');

    models{i} = m;
    fprintf('i=%d;\n',i);
end