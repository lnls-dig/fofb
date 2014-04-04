function simout = simfofb_sirius

% Matrices
nbpm = 5;
ncorr = 4;
ndist = 2;
param.RespMat = floor(10*rand(nbpm, ncorr));
param.CorrMat = calcpseudoinv(param.RespMat);
param.DistMat = floor(10*rand(nbpm, ndist));

% BPM (FOFB sensor)
bpm_cic_decfactor = 1250;
bpm_cic_nstages = 3;
bpm_cic_ndelays = 1;
[param.bpm.num, param.bpm.den] = buildcic(bpm_cic_decfactor, bpm_cic_nstages, bpm_cic_ndelays);
param.bpm.Ts = 10e-6;

% DCCT (Power supply sensor)
param.dcct.num = 1;
param.dcct.den = 1;

% Power supply controller
param.psctrl.num = 1;
param.psctrl.den = 1;
param.psctrl.Ts = 10e-6;
param.psctrl.manual = 0;   % (0 = auto; 1 = manual)

% Power supply output filter
param.psfilter.num = 1;
param.psfilter.den = 1;

% Corrector magnet
corrmagnet_L = 20e-3;
corrmagnet_R = 1;
param.corrmagnet.num = 1;
param.corrmagnet.den = [corrmagnet_L/corrmagnet_R 1];

% Vacuum chamber
vac_radius = 12e-3;
vac_thickness = 0.1e-3;
vac_conductivity = 5.85e7;
mu0 = 4*pi*1e-7;
tau = 0.5*mu0*vac_conductivity*vac_radius*vac_thickness;
param.vacchamb.num = 1;
param.vacchamb.den = [tau 1];

% FOFB controller
[a,b,c,d] = tf2ss(0.0001*[1 0], [1 -1]);
[param.fofbctrl.a, b, param.fofbctrl.c, d] = repss(a,b,c,d,ncorr);
param.fofbctrl.b = b*param.CorrMat;
param.fofbctrl.d = d*param.CorrMat;
param.fofbctrl.Ts = 10e-6;
param.fofbctrl.sensdelay = 10e-6;
param.fofbctrl.actdelay = 10e-6;

% Build simulation vectors
t = (0:1e-6:1)';
npts = length(t);
param.simvectors.t = t;
param.simvectors.dist = 0.001*[zeros(200000,ndist); ones(npts-200000,ndist)];
param.simvectors.bpm_noise = 140e-9*randn(npts,nbpm);
param.simvectors.dcct_noise = 3e-4*randn(npts,ncorr);
param.simvectors.psoutput_noise = 3e-4*randn(npts,ncorr);

simout = simfofb(param);