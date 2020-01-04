function simout = simofbuvx

% Load LNLS-UVX storage ring matrices (Horizontal Plane)
try
    load RespMatH_LNLS
    load RespDistH_LNLS
catch
    error('Matrix files ''RespMatH_LNLS.mat'' and ''RespDistH_LNLS.mat'' must be in the path.');
end

param.RespMat = RespMatH_LNLS;
param.CorrMat = pseudoinv(RespMatH_LNLS);
param.DistMat = RespDistH_LNLS;

% Dimensions
nbpm = size(param.RespMat, 1);
ncorr = size(param.RespMat, 2);
ndist = size(param.DistMat, 2);

% BPM (FOFB sensor) SISO model
[param.bpm.num, param.bpm.den] = butter(2, 1/3);
param.bpm.Ts = 330e-6;

% Corrector magnet SISO model
param.corrmagnet.num = 1;
param.corrmagnet.den = [700e-6 1];

% Vacuum chamber SISO model
param.vacchamb.num = 1;
param.vacchamb.den = [1/2/pi/1250 1];

% FOFB controller SISO model
[a,b,c,d] = tf2ss(0.05, [1 -1]);
[param.fofbctrl.a, b, param.fofbctrl.c, d] = repss(a,b,c,d,ncorr);
param.fofbctrl.b = b*param.CorrMat;
param.fofbctrl.d = d*param.CorrMat;
param.fofbctrl.Ts = 330e-6;
param.fofbctrl.sensdelay = 330e-6;
param.fofbctrl.actdelay = 330e-6;

% DCCT SISO model (Power supply sensor) (feed-through)
param.dcct.num = 1;
param.dcct.den = 1;

% Power supply controller SISO model (feed-through)
param.psctrl.num = 1;
param.psctrl.den = 1;
param.psctrl.Ts = 330e-6;
param.psctrl.manual = 1;

% Power supply output filter SISO model (feed-through)
param.psfilter.num = 1;
param.psfilter.den = 1;

% Build simulation vectors
t = (0:100e-6:1);
npts = length(t);
param.simvectors.t = t(:);
param.simvectors.dist = 1e-3*[zeros(2000,ndist); ones(npts-2000,ndist)];
param.simvectors.bpm_noise = 10e-6*randn(npts,nbpm);
param.simvectors.dcct_noise = 10e-6*randn(npts,ncorr);
param.simvectors.psoutput_noise = 10e-6*randn(npts,ncorr);

% Simulate FOFB MIMO model
simout = simofb(param);