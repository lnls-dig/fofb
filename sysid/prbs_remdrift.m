function [y, nseg_discard] = prbs_remdrift(y, prbslen)

fir_remove_drift = dfilt.dffir(ones(1, prbslen)/prbslen);
drift_y = filter(fir_remove_drift, y);

dly = grpdelay(fir_remove_drift);
dly = dly(1);
dly = ceil(dly);

drift_y = [drift_y(1+dly:end, :); zeros(dly, size(y,2))];

y = y-drift_y;
nseg_discard = ceil(dly/prbslen);