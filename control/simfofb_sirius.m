function simout = simfofb_sirius

load ../matlab.mat

M = M(:,ceil((end+1)/2):end,3);
Md = Mdq(:,:,3);


% Matrices
nbpm = size(M,1);
ncorr = size(M,2);
ndist = size(Md,2);
param.RespMat = M;
param.CorrMat = pinv(param.RespMat);
param.DistMat = Md;

% BPM (FOFB sensor)
bpm_cic_decfactor = 1015;
bpm_cic_nstages = 2;
bpm_cic_ndelays = 3;
[param.bpm.num, param.bpm.den] = buildcic(bpm_cic_decfactor, bpm_cic_nstages, bpm_cic_ndelays);
param.bpm.Ts = 10e-6;

% DCCT (Power supply sensor)
param.dcct.num = 1;
param.dcct.den = 1;

% Power supply controller
param.psctrl.num = 1;
param.psctrl.den = 1;
param.psctrl.Ts = 10e-6;
param.psctrl.manual = true;

% Power supply output filter
param.psfilter.num = 1;
param.psfilter.den = 1;

% Corrector magnet
corrmagnet_L = 20e-3;
corrmagnet_R = 1;
[param.corrmagnet.num, param.corrmagnet.den] = build1order(corrmagnet_L/corrmagnet_R);

% Vacuum chamber
vac_radius = 12e-3;
vac_thickness = 0.1e-3;
vac_conductivity = 5.85e7;
mu0 = 4*pi*1e-7;
[param.vacchamb.num, param.vacchamb.den] = build1order(0.5*mu0*vac_conductivity*vac_radius*vac_thickness);

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