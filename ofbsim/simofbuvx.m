function simout = simofbuvx

% Response matrices
load matrix_lnls1.mat
respmat = Mh_meas;
distmat = Mdisth_teor;

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
[a,b,c,d] = tf2ss(num, den);
[a,b,c,d] = repss(a,b,c,d,ncorrfofb);
param.fofbctrl = ss(a, b*corrmat, c, d*corrmat, -1);

% BPM model
[num, den] = buildcic(32, 2, 1);
param.bpm = tf(num, den, -1);

% Corrector magnet current model
corrmagnet_L = 0.4e-3;
corrmagnet_R = 1;
[num, den] = build1order(corrmagnet_L/corrmagnet_R);
param.corrmagnet = tf(num, den);

% Corrector magnet current/magnetic field gain
param.corrmagnetgain = ones(ncorrfofb,1);

% Vacuum chamber model
[num, den] = build1order(1/1.25e3/2/pi);
param.vacchamb = tf(num, den);

% Communication network delays
param.netdelay.bpm = Ts;
param.netdelay.corrfofb = Ts;

% Sensed orbit points (orbit measured by BPMs)
param.bpmordering = 1:nbpm;

% Simulation vectors
t = (0:10e-6:Ts*1000)';
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