function simout = simfofb_sirius

load ../matlab.mat

respmat = M(:,ceil((end+1)/2):end,3);
distmat = Mdq(:,:,3);

% Matrices
nbpm = size(respmat,1);
ncorr = size(respmat,2);
ndist = size(distmat,2);
corrmat = pinv(respmat);

% Beam dynamics MIMO model
param.beam.sys = ss([],[],[],respmat);

% Beam disturbance MIMO model
param.distbeam.sys = ss([],[],[],distmat);

% BPM (FOFB sensor)
bpm_cic_decfactor = 1015;
bpm_cic_nstages = 2;
bpm_cic_ndelays = 1;
[param.bpm.num, param.bpm.den] = buildcic(bpm_cic_decfactor, bpm_cic_nstages, bpm_cic_ndelays);
param.bpm.Ts = 10e-6;

% Power supply output filter
param.corrsofbpsfilter.num = 1;
param.corrsofbpsfilter.den = 1;

% Corrector magnet current
corrmagnet_L = 20e-3;
corrmagnet_R = 1;
[param.corrmagnet.num, param.corrmagnet.den] = build1order(corrmagnet_L/corrmagnet_R);

% Corrector magnet current/magnetic field gain
param.corrmagnetgain = ones(ncorr,1);

% Vacuum chamber
vac_radius = 12e-3;
vac_thickness = 0.1e-3;
vac_conductivity = 5.85e7;
mu0 = 4*pi*1e-7;
[param.vacchamb.num, param.vacchamb.den] = build1order(0.5*mu0*vac_conductivity*vac_radius*vac_thickness);

% FOFB controller
num = 0.001*[1 0];
den = [1 -1];
[num,den] = eqtflength(num, den);
[a,b,c,d] = tf2ss(num, den);
[a,b,c,d] = repss(a,b,c,d,ncorr);
param.fofbctrl.sys = ss(a, b*corrmat, c, d*corrmat, 10e-6);

% SOFB controller
param.sofbctrl.sys = ss([], [], [], ones(1,nbpm), 1);

% Communication network delays
param.netdelay.bpm = 11e-6;
param.netdelay.corrfofb = 12e-6;
param.netdelay.corrsofb = 103e-6;


%%%%%%%%%%%%
ncorrsofb = 1;
ncorrfofb = size(param.fofbctrl.sys.d, 1);
param.corrselectsofb = 1;
param.ncorr = 1;
param.corrordering = 1+(1:ncorr);


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