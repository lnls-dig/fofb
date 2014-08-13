function [param, fit] = idgbbpm(datae, datav, order_bpm, order_corr, abpm, bbpm, acorr, bcorr, M, bpmungain, bpmhvsymmetry)

[~, nbpm, ncorr] = size(datae);

param = [idab2param(abpm, bbpm, order_bpm, bpmungain) idab2param(acorr, bcorr, order_corr)];

auxvars = [nbpm, ncorr, order_bpm, order_corr, [0 0], bpmungain, 0, bpmhvsymmetry, M(:)'];
fofbmodel = idgrey('idfofbmodel', param, 'd', auxvars, 'InitialState', 'Zero', 'OutputName', datav.OutputName, 'InputName', datav.InputName);
fofbmodel.Ts = datae.Ts;
fofbmodel.Algorithm.Focus = 'Simulation';
fofbmodel.Algorithm.Display = 'On';

fofbmodele = pem(detrend(datae,0), fofbmodel);

datav_ = detrend(datav,0);
[~,fit] = compare(fofbmodele, datav_, Inf);
fit = squeeze(fit);

param = fofbmodele.ParameterVector;