function simout = simfofb(param)
%SIMFOFB Prepare structure arrays for a Simulink-based simualtion of a Fast
%   Orbit Feedback loop (FOFB) and call the Simulink model 'fofb.mdl'.
%
%   simout = simfofb(param)
%   
%   Inputs
%       PARAM: structure array containing all necessary information for
%           simulating the FOFB model, comprising FOFB controller,
%           SISO models of all dynamics systems, orbit matrices and
%           simulation vectors. The structure of PARAM is as follows:
%               param
%                 .RespMat: Orbit Response Matrix (orbit vs. orbit corrector)
%                 .CorrMat: Orbit Correction Matrix (orbit corrector vs. orbit)
%                 .DistMat: Orbit Disturbances Matrix (orbit vs. source of orbit disturbance)
%                 .fofbctrl: digital state-space representation of the FOFB controller (MIMO model)
%                   .a
%                   .b
%                   .c
%                   .d
%                   .Ts
%                 .bpm: digital transfer function of the Beam Position Monitor (BPM) - FOFB sensor (SISO model)
%                   .num
%                   .den
%                   .Ts
%                 .psctrl: digital transfer fucntion of the orbit corrector power supply (SISO model)
%                   .num
%                   .den
%                   .Ts
%                 .dcct: analog transfer function of the DC Current Transformer - power supply sensor (SISO model)
%                   .num
%                   .den
%                 .psfilter: analog transfer function of the orbit corrector power supply's output filter (SISO model)
%                   .num
%                   .den
%                 .corrmagnet: analog transfer function of the orbit corrector's magnet (SISO model)
%                   .num
%                   .den
%                 .vacchamb: analog transfer function of the vacuum chamber eddy currents effect (SISO model)
%                   .num
%                   .den
%                 .simvectors: simulation vectors
%                   .t: time vector (n x 1 or 1 x n array)
%                   .dist: disturbance vectors (n x <number of source of disturbances> matrix)
%                   .bpm_noise: BPM noise vectors (n x <number of BPMs> matrix)
%                   .dcct_noise: DCCT noise vectors (n x <number of orbit correctors> matrix)
%                   .psoutput_noise: power supply output noise vectors (n x <number of orbit correctors> matrix)
%
%   Outputs
%       SIMOUT: structure array containing simulation results. The
%           structure of SIMOUT is as follows:
%               simout
%                   .t: simulation time
%                   .beam: real beam position
%                   .control_fofb: FOFB controller's control signal output
%                   .sensor_fofb: FOFB controller's sensor signal input
%                   .control_ps: power supply controller's control signal output
%                   .sensor_ps: power supply controller's sensor signal input

RespMat = param.RespMat;
CorrMat = param.CorrMat;
DistMat = param.DistMat;
bpm = param.bpm;
dcct = param.dcct;
vacchamb = param.vacchamb;
fofbctrl = param.fofbctrl;
psctrl = param.psctrl;
psfilter = param.psfilter;
corrmagnet = param.corrmagnet;
simvectors = param.simvectors;

% Check if all matrix dimensions are correct
if size(RespMat, 1) ~= size(CorrMat, 2)
    error('CorrMat must have the same number of columns as the number of rows of RespMat (number of BPMs).');
elseif size(RespMat, 2) ~= size(CorrMat, 1)
    error('CorrMat must have the same number of rows as the number of columns of RespMat (number of orbit correctors).');
elseif size(RespMat, 1) ~= size(DistMat, 1)
    error('DistMat must have the same number of rows as the number of rows of RespMat (number of BPMs).');
end

% Extract number of BPMs, orbit correctors and sources of disturbance
nbpm = size(RespMat, 1);
ncorr = size(RespMat, 2);
ndist = size(DistMat, 2);

% Check if all simulation vector dimensions are correct
if size(simvectors.dist, 2) ~= ndist
    error('simvectors.dist must have the same number of columns as the number of columns of DistMat (number of source of disturbances).');
elseif size(simvectors.bpm_noise, 2) ~= nbpm
    error('simvectors.bpm_noise must have the same number of columns as the number of rows of RespMat (number of BPMs).');
elseif size(simvectors.dcct_noise, 2) ~= ncorr
    error('simvectors.dcct_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
elseif size(simvectors.psoutput_noise, 2) ~= ncorr
    error('simvectors.psoutput_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
elseif size(simvectors.dist, 1) ~= length(simvectors.t)
    error('simvectors.dist must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.bpm_noise, 1) ~= length(simvectors.t)
    error('simvectors.bpm_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.dcct_noise, 1) ~= length(simvectors.t)
    error('simvectors.dcct_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.psoutput_noise, 1) ~= length(simvectors.t)
    error('simvectors.psoutput_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
end

% Beam dynamics MIMO model
beam.a = [];
beam.b = [];
beam.c = [];
beam.d = RespMat;

% Beam disturbance MIMO model
distbeam.a = [];
distbeam.b = [];
distbeam.c = [];
distbeam.d = DistMat;

% BPM  MIMO model (FOFB sensors)
[a,b,c,d] = tf2ss(bpm.num, bpm.den);
[bpm.a, bpm.b, bpm.c, bpm.d] = repss(a,b,c,d,nbpm);

% DCCT MIMO model (Power supplies sensors)
[a,b,c,d] = tf2ss(dcct.num, dcct.den);
[dcct.a, dcct.b, dcct.c, dcct.d] = repss(a,b,c,d,ncorr);

% Power supply controller MIMO model
[a,b,c,d] = tf2ss(psctrl.num, psctrl.den);
[psctrl.a, psctrl.b, psctrl.c, psctrl.d] = repss(a,b,c,d,ncorr);
psctrl.Ts = psctrl.Ts;

% Power supply output filter MIMO model
[a,b,c,d] = tf2ss(psfilter.num, psfilter.den);
[psfilter.a, psfilter.b, psfilter.c, psfilter.d] = repss(a,b,c,d,ncorr);

% Corrector magnet MIMO model
[a,b,c,d] = tf2ss(corrmagnet.num, corrmagnet.den);
[corrmagnet.a, corrmagnet.b, corrmagnet.c, corrmagnet.d] = repss(a,b,c,d,ncorr);

% Vacuum chamber MIMO model
[a,b,c,d] = tf2ss(vacchamb.num, vacchamb.den);
[vacchamb.a, vacchamb.b, vacchamb.c, vacchamb.d] = repss(a,b,c,d,ncorr);

% Assign simulation parameters to workspace so that Simulink can use it
simparam = struct(...
    'beam', beam, ...
    'distbeam', distbeam, ...
    'fofbctrl', fofbctrl, ...
    'bpm', bpm, ...
    'dcct', dcct, ...
    'psctrl', psctrl, ...
    'psfilter', psfilter, ...
    'corrmagnet', corrmagnet, ...
    'vacchamb', vacchamb, ...
    'simvectors', simvectors);

assignin('base', 'simparam', simparam);

% Run Simulink simulation
simout_simulink = sim('fofb', 'StopTime', num2str(simvectors.t(end)));

% Extract simulation results
simout.t = get(simout_simulink, 'tout');
yout = get(simout_simulink, 'yout');
simout.beam = yout(:,1:nbpm);
simout.control_fofb = yout(:,nbpm+(1:ncorr));
simout.control_ps = yout(:,ncorr+nbpm+(1:ncorr));
simout.sensor_fofb = yout(:,2*ncorr+nbpm+(1:nbpm));
simout.sensor_ps = yout(:,2*ncorr+2*nbpm+(1:ncorr));