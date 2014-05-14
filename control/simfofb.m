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

respmat = param.beam.sys.d;
distmat = param.distbeam.sys.d;


beam = param.beam;
distbeam = param.distbeam;
bpm = param.bpm;
sofbctrl = param.sofbctrl;
fofbctrl = param.fofbctrl;

corrsofbpsfilter = param.corrsofbpsfilter;
corrmagnet = param.corrmagnet;
corrmagnetgain = param.corrmagnetgain;
vacchamb = param.vacchamb;
netdelay = param.netdelay;
simvectors = param.simvectors;

% Check if all matrix dimensions are correct
if size(respmat, 1) ~= size(distmat, 1)
    error('distmat must have the same number of rows as tpsctrl.Ts = psctrl.Ts;he number of rows of respmat (number of beam position readings).');
end

% Extract number of BPMs, orbit correctors and sources of disturbance
nbpm = size(respmat, 1);
ncorr = size(respmat, 2);
ncorrfofb = size(fofbctrl.sys.d, 1);
ncorrsofb = size(sofbctrl.sys.d, 1);
ndist = size(distmat, 2);

% check ncorr == ncorrfofb+ncorrsofb

% Check if all simulation vector dimensions are correct
if size(simvectors.orbit_distortion, 2) ~= ndist
    error('simvectors.orbit_distortion must have the same number of columns as the number of columns of DistMat (number of source of disturbances).');
elseif size(simvectors.bpm_noise, 2) ~= nbpm
    error('simvectors.bpm_noise must have the same number of columns as the number of rows of RespMat (number of BPMs).');
elseif size(simvectors.corrsofb_current_noise, 2) ~= ncorrsofb
    error('simvectors.corrsofb_current_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
elseif size(simvectors.corrsofb_power_noise, 2) ~= ncorrsofb
    error('simvectors.corrsofb_power_noise must have the same number of columns as the number of columns of RespMat (number of orbit correctors).');
elseif size(simvectors.orbit_distortion, 1) ~= length(simvectors.t)
    error('simvectors.orbit_distortion must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.bpm_noise, 1) ~= length(simvectors.t)
    error('simvectors.bpm_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.corrsofb_current_noise, 1) ~= length(simvectors.t)
    error('simvectors.corrsofb_current_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
elseif size(simvectors.corrsofb_power_noise, 1) ~= length(simvectors.t)
    error('simvectors.corrsofb_power_noise must have the same number of rows as the number of elements of simvectors.t (number of simulation samples).');
end

% Communication network delays
fracnetdelay.bpm      = repmat(rem(netdelay.bpm, fofbctrl.sys.Ts), 1, nbpm);
fracnetdelay.corrfofb = repmat(rem(netdelay.corrfofb, fofbctrl.sys.Ts), 1, ncorrfofb);
fracnetdelay.corrsofb = repmat(rem(netdelay.corrsofb, sofbctrl.sys.Ts), 1, ncorrsofb);

% BPM MIMO model
[num,den] = eqtflength(bpm.num, bpm.den);
[a,b,c,d] = tf2ss(num, den);
[a,b,c,d] = repss(a,b,c,d,nbpm);
bpm.sys = ss(a,b,c,d, fofbctrl.sys.Ts);
delay = floor(fracnetdelay.bpm/fofbctrl.sys.Ts);
%if delay > 0
    bpm.sys.OutputDelay = 1+0*delay;
%end

% Slow corrector power supply output filter MIMO model
[a,b,c,d] = tf2ss(corrsofbpsfilter.num, corrsofbpsfilter.den);
[a,b,c,d] = repss(a,b,c,d,ncorrsofb);
corrsofbpsfilter.sys = ss(a,b,c,d);

% Corrector magnets MIMO model
[a,b,c,d] = tf2ss(corrmagnet.num, corrmagnet.den);
[a,b,c,d] = repss(a,b,c,d,ncorr);
corrmagnet.sys = ss(a,b,c,d);

% Vacuum chamber MIMO model
[a,b,c,d] = tf2ss(vacchamb.num, vacchamb.den);
[a,b,c,d] = repss(a,b,c,d,ncorr);
vacchamb.sys = ss(a,b,c,d);


delay = floor(fracnetdelay.corrfofb/fofbctrl.sys.Ts);
if delay > 0
    fofbctrl.sys.OutputDelay = delay;
end

delay = floor(fracnetdelay.corrsofb/sofbctrl.sys.Ts);
if delay > 0
    sofbctrl.sys.OutputDelay = delay;
end

% Quantization
% quantization.bpm = repmat(quantization.bpm, 1, nbpm);
% quantization.corrfofb = repmat(quantization.corrfofb, 1, ncorrfofb);
% quantization.corrsofb = repmat(quantization.corrsofb, 1, ncorrsofb);
%assignin('base', 'quantization', quantization);

modelname = 'ofb';
open_system(modelname,'loadonly');
if true
    set_param([modelname '/Slow Orbit Correctors'' Power Supplies'], 'ModelFile', 'corrsofb_ps_voltage.mdl');
elseif strcmpi(sofb_pstype, 'x')
    pssens = param.psdcct;
    psctrl = param.psctrl;

    % Slow corrector power supply sensors MIMO model
    [a,b,c,d] = tf2ss(pssens.num, pssens.den);
    [pssens.a, pssens.b, pssens.c, pssens.d] = repss(a,b,c,d,ncorr);
    
    % Slow corrector power supply controller MIMO model
    [a,b,c,d] = tf2ss(psctrl.num, psctrl.den);
    [psctrl.a, psctrl.b, psctrl.c, psctrl.d] = repss(a,b,c,d,ncorr);
    
    assignin('base', 'pssens', pssens);
    assignin('base', 'psctrl', psctrl);
elseif strcmpi(sofb_pstype, 'y')
end
save_system(modelname);
%close_system(modelname);

% Assign simulation parameters to workspace so that Simulink can use it
assignin('base', 'beam', beam);
assignin('base', 'distbeam', distbeam);
assignin('base', 'sofbctrl', sofbctrl);
assignin('base', 'fofbctrl', fofbctrl);
assignin('base', 'bpm', bpm);
assignin('base', 'corrsofbpsfilter', corrsofbpsfilter);
assignin('base', 'corrmagnetgain', corrmagnetgain);
assignin('base', 'corrmagnet', corrmagnet);
assignin('base', 'vacchamb', vacchamb);
assignin('base', 'fracnetdelay', fracnetdelay);
assignin('base', 'nbpm', nbpm);
assignin('base', 'ndist', ndist);
assignin('base', 'ncorr', ncorr);
assignin('base', 'ncorrsofb', ncorrsofb);
assignin('base', 'ncorrfofb', ncorrfofb);
assignin('base', 'simvectors', simvectors);
assignin('base', 'corrselectsofb', param.corrselectsofb);
assignin('base', 'corrordering', param.corrordering);

% Run Simulink simulation
simout_simulink = sim(modelname, 'StopTime', num2str(simvectors.t(end)));

% Extract simulation results
r = get(simout_simulink, 'yout');
simout.t = r.time;
simout.beam_orbit = r.signals(1).values;
simout.control_sofb = r.signals(2).values;
simout.control_fofb = r.signals(3).values;
simout.corr_voltages = r.signals(4).values;
simout.corr_currents = r.signals(5).values;