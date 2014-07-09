function simout = simofbuvx

% Response matrices
load matrix_lnls1.mat
respmat =  Mv_meas;
distmat =  ones(25,1); Mdisth_teor;

% Sampling period
Ts = 320e-6;
param.Ts = Ts;

% Dimensions
norbit = size(respmat,1);
nbpm = norbit;
ncorrfofb = size(respmat,2);
ndist = size(distmat,2);

% Beam orbit model
param.beam = ss([],[],[],respmat);

% Beam disturbance model
param.distbeam = ss([],[],[],distmat);

% FOFB controller model
corrmat = pinv(respmat);
num = 0.05*[1 0];
den = [1 -1];
%num = [1.12007392895013 -1.48080307715986 0.462362520083246];
%den = [1 -1 0];
[a,b,c,d] = tf2ss(num, den);
[a,b,c,d] = repss(a,b,c,d,ncorrfofb);
param.fofbctrl = ss(a, b*corrmat, c, d*corrmat, -1);

% BPM model
[num, den] = buildcic(32, 2, 1);
param.bpm = tf(num, den, -1);

% Corrector magnet current model
param.corrmagnet_current = corrmagnetv(1,1);
param.corrmagnet_current.inputdelay = 0;

% Corrector magnet current/magnetic field gain
param.corrmagnetgain = ones(ncorrfofb,1);

% % Corrector magnet core
% corrmagnet_core = zpk(-1094,[-991.1 -9005],8162.1842);
% param.corrmagnet_core = corrmagnet_core;

%8162.1842 (s+1094)
%------------------
%(s+991.1) (s+9005)

% Vacuum chamber model
[num, den] = build1order(1/2.5e3/2/pi);
vacchamb = tf(num, den);
param.vacchamb = vacchamb;

% Communication network delays
param.netdelay.bpm = Ts;
param.netdelay.corrfofb = Ts;

% Sensed orbit points (orbit measured by BPMs)
param.bpmordering = 1:nbpm;

% Simulation vectors
t = (0:32e-6:Ts*10000)';
npts = length(t);
param.simvectors.t = t;
param.simvectors.reference_orbit = zeros(npts,nbpm);
param.simvectors.orbit_distortion = [zeros(2000,ndist); ones(npts-2000,ndist)];
param.simvectors.fofb_setpoints = zeros(npts, ncorrfofb);
param.simvectors.bpm_noise = 0*140e-9*randn(npts,nbpm);
param.simvectors.sofb_setpoints = zeros(npts, 0);
param.simvectors.corrsofb_current_noise = 3e-4*randn(npts, 0);
param.simvectors.corrsofb_power_noise = 3e-4*randn(npts, 0);

simout = simofb(param);