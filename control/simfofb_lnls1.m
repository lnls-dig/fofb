function simout = simfofb_lnls1

load matrix_lnls1.mat

plane = 'h';

if strcmpi(plane, 'h')
    respmat = Mh_meas;
    distmat = Mdisth_teor;
else
    respmat = Mv_meas;
    distmat = Mdistv_teor;
end

distmat = ones(size(respmat,1),1);

% Matrices
norbit = size(respmat,1);
ncorr = size(respmat,2);
ndist = size(distmat,2);
corrmat = pinv(respmat);

nbpm = norbit;

% Beam dynamics MIMO model
param.beam.sys = ss([],[],[],respmat);

% Beam disturbance MIMO model
param.distbeam.sys = ss([],[],[],distmat);

% BPM (FOFB sensor)
bpm_cic_decfactor = 1015;
bpm_cic_nstages = 2;
bpm_cic_ndelays = 1;
[param.bpm.num, param.bpm.den] = buildcic(bpm_cic_decfactor, bpm_cic_nstages, bpm_cic_ndelays);

% Power supply output filter
param.corrsofbpsfilter.num = 1;
param.corrsofbpsfilter.den = 1;

% Corrector magnet current
corrmagnet_L = 0.4e-3;
corrmagnet_R = 1;
[param.corrmagnet.num, param.corrmagnet.den] = build1order(corrmagnet_L/corrmagnet_R);

% Corrector magnet current/magnetic field gain
param.corrmagnetgain = ones(ncorr,1);

% Vacuum chamber
[param.vacchamb.num, param.vacchamb.den] = build1order(1/1.25e3/2/pi);

% FOFB controller
num = 0.005*[1 0];
den = [1 -1];
[num,den] = eqtflength(num, den);
[a,b,c,d] = tf2ss(num, den);
[a,b,c,d] = repss(a,b,c,d,ncorr);
param.fofbctrl.sys = ss(a, b*corrmat, c, d*corrmat, 320e-6);

% SOFB controller
param.sofbctrl.sys = ss([], [], [], ones(1,nbpm), 1);

% Communication network delays
param.netdelay.bpm = 320e-6;
param.netdelay.corrfofb = 320e-6;
param.netdelay.corrsofb = 0e-6;


%%%%%%%%%%%%
ncorrsofb = 1;
ncorrfofb = size(param.fofbctrl.sys.d, 1);
param.corrsofbselect = 1;
param.ncorr = 1;
param.corrordering = 1+(1:ncorr);
param.bpmselect = 1:nbpm;

% Build simulation vectors
t = (0:1e-6:0.1)';
npts = length(t);
param.simvectors.t = t;
param.simvectors.reference_orbit = zeros(npts,nbpm);
param.simvectors.orbit_distortion = 0.001*[zeros(2000,ndist); ones(npts-2000,ndist)];
param.simvectors.sofb_setpoints = zeros(npts, ncorrsofb);
param.simvectors.fofb_setpoints = zeros(npts, ncorrfofb);
param.simvectors.bpm_noise = 140e-9*randn(npts,nbpm);
param.simvectors.corrsofb_current_noise = 3e-4*randn(npts,ncorrsofb);
param.simvectors.corrsofb_power_noise = 3e-4*randn(npts,ncorrsofb);

simout = simfofb(param);