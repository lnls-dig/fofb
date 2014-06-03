function simout = simofb(param)
%SIMOFB Prepare structure arrays for a Simulink-based simualtion of a Fast
%   Orbit Feedback loop (FOFB) and call the Simulink model 'fofb.mdl'.
%
%   simout = simofb(param)
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
%                 .bpm: digital transfer function of the Beam Position Monitor (BPM) - FOFB sensor (SISO model)
%                 .psctrl: digital transfer fucntion of the orbit corrector power supply (SISO model)
%                 .pssens: analog transfer function of the DC Current Transformer - power supply sensor (SISO model)
%                 .psfilter: analog transfer function of the orbit corrector power supply's output filter (SISO model)
%                 .corrmagnet: analog transfer function of the orbit corrector's magnet (SISO model)
%                 .vacchamb: analog transfer function of the vacuum chamber eddy currents effect (SISO model)
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

modelname = 'ofb';

respmat = param.beam.d;
distmat = param.distbeam.d;

beam = param.beam;
distbeam = param.distbeam;
bpm = param.bpm;
fofbctrl = param.fofbctrl;
corrmagnet = param.corrmagnet;
corrmagnetgain = param.corrmagnetgain;
vacchamb = param.vacchamb;
netdelay = param.netdelay;

simvectors = param.simvectors;
Ts = param.Ts;

% Check if all matrix dimensions are correct
if size(respmat, 1) ~= size(distmat, 1)
    error('distmat must have the same number of rows as tpsctrl.Ts = psctrl.Ts;he number of rows of respmat (number of beam position readings).');
end

% Dimensions
norbit = size(respmat, 1);
nbpm = size(fofbctrl.d, 2);
ncorrfofb = size(fofbctrl.d, 1);
ndist = size(distmat, 2);
npts = size(param.simvectors.t, 1);

open_system(modelname,'loadonly');
if ~isfield(param, 'sofbctrl') || isempty(param.sofbctrl) || ~isdt(param.sofbctrl)
    ncorrsofb = 1;
    s2ffactor = 1;
    simvectors.sofb_setpoints = zeros(npts,1);
    simvectors.corrsofb_current_noise = zeros(npts,1);
    simvectors.corrsofb_power_noise = zeros(npts,1);
    set_param([modelname '/SOFB'], 'ModelFile', 'sofb_bypass.mdl');
    ncorr = ncorrfofb;
else
    sofbctrl = param.sofbctrl;
    ncorrsofb = size(sofbctrl.d, 1);
    s2ffactor = param.s2ffactor;
    assignin('base', 'sofbctrl', sofbctrl);
    set_param([modelname '/SOFB'], 'ModelFile', 'sofb.mdl');
    ncorr = ncorrfofb + ncorrsofb;
end
save_system(modelname);

% Check if all simulation vector dimensions are correct
% TODO: need to finalize!
% check ncorr == ncorrfofb+ncorrsofb

% check discrete models with isdt and sys.Ts == -1

% if size(simvectors.orbit_distortion, 2) ~= ndist
%     error('simvectors.orbit_distortion must have the same number of columns as the number of columns of DistMat (number of source of disturbances).');
% elseif size(simvectors.bpm_noise, 2) ~= nbpm
%     error('simvectors.bpm_noise must have the same number of columns as the number of rows of RespMat (number of BPMs).');
% elseif size(simvectors.corrsofb_current_noise, 2) ~= ncorrsofb
%     error('simvectors.corrsofb_current_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
% elseif size(simvectors.corrsofb_power_noise, 2) ~= ncorrsofb
%     error('simvectors.corrsofb_power_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
% elseif size(simvectors.orbit_distortion, 1) ~= length(simvectors.t)
%     error('simvectors.orbit_distortion must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
% elseif size(simvectors.bpm_noise, 1) ~= length(simvectors.t)
%     error('simvectors.bpm_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
% elseif size(simvectors.corrsofb_current_noise, 1) ~= length(simvectors.t)
%     error('simvectors.corrsofb_current_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
% elseif size(simvectors.corrsofb_power_noise, 1) ~= length(simvectors.t)
%     error('simvectors.corrsofb_power_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
% end

% BPM MIMO model
if length(bpm) == 1
    bpm = siso2mimo(bpm, nbpm);
end

% Corrector magnets MIMO model
if length(corrmagnet) == 1
    corrmagnet = siso2mimo(corrmagnet, ncorr);
end

% Vacuum chamber MIMO model
if length(vacchamb) == 1
    vacchamb = siso2mimo(vacchamb, ncorr);
end

% Delays (integer part - embed inside LTI systems)
delay.bpm = ss([],[],[],eye(nbpm),-1,'OutputDelay', floor(netdelay.bpm/Ts));
delay.corrfofb = ss([],[],[],eye(ncorrfofb),-1,'OutputDelay', floor(netdelay.corrfofb/Ts));

% Delays (fractional part - use Simulink's transport delay block)
fracdelay.bpm      = rem(netdelay.bpm, Ts);
fracdelay.corrfofb = rem(netdelay.corrfofb, Ts);

% Quantization
% quantization.bpm = repmat(quantization.bpm, 1, nbpm);
% quantization.corrfofb = repmat(quantization.corrfofb, 1, ncorrfofb);
% quantization.corrsofb = repmat(quantization.corrsofb, 1, ncorrsofb);
%assignin('base', 'quantization', quantization);

% Assign simulation parameters to workspace so that Simulink can use it
assignin('base', 'beam', beam);
assignin('base', 'distbeam', distbeam);
assignin('base', 'fofbctrl', fofbctrl);
assignin('base', 'bpm', bpm);
assignin('base', 'corrmagnetgain', corrmagnetgain);
assignin('base', 'corrmagnet', corrmagnet);
assignin('base', 'vacchamb', vacchamb);
assignin('base', 'delay', delay);
assignin('base', 'fracdelay', fracdelay);
assignin('base', 'norbit', norbit);
assignin('base', 'ndist', ndist);
assignin('base', 'nbpm', nbpm);
assignin('base', 'ncorrsofb', ncorrsofb);
assignin('base', 'ncorrfofb', ncorrfofb);
assignin('base', 'ncorrmagnet', ncorr);
assignin('base', 'simvectors', simvectors);
assignin('base', 'bpmordering', param.bpmordering);
assignin('base', 's2ffactor', s2ffactor);
assignin('base', 'Ts', Ts);

% Run simulation model (Simulink)
simout_simulink = sim(modelname, 'StopTime', num2str(simvectors.t(end)));

% Extract simulation results
r = get(simout_simulink, 'yout');
simout.t = r.time;
simout.beam_orbit = r.signals(1).values;
simout.beam_orbit_bpm = r.signals(2).values;
simout.control_fofb = r.signals(3).values;
simout.control_sofb = r.signals(4).values;
simout.corr_voltages = r.signals(5).values;
simout.corr_currents = r.signals(6).values;