function [bpm, corr] = fofb_sisomodel(config)

% Controller sample time
Ts = 1/config.ctrl_sampling_rate;

% BPM electronics (CIC decimation)
cic_nstages = config.bpm_cic_nstages;
cic_ndelays = config.bpm_cic_ndelays;
cic_decimation_factor = config.bpm_cic_decimation_factor;
cic_coeff = tf(mfilt.cicdecim(cic_decimation_factor, cic_ndelays, cic_nstages));
decimated_cic_coeff = interp1(1:length(cic_coeff), cic_coeff, linspace(1,length(cic_coeff), cic_nstages*cic_ndelays + 1));
decimated_cic_coeff((ceil(length(decimated_cic_coeff)/2)+1):end) = decimated_cic_coeff((floor(length(decimated_cic_coeff)/2)):-1:1);
decimated_cic_coeff = decimated_cic_coeff/sum(decimated_cic_coeff);
bpm = tf(decimated_cic_coeff,[1 zeros(1,length(decimated_cic_coeff)-1)],Ts);

% Delays on data distribution + FOFB control algorithm processing
latency = config.communication_latency + Ts;

% Actuator (power supply + orbit corrector magnet)
if strcmpi(config.actuator, 'voltage')
    % L-R model for orbit corrector magnet (first-order model)
    cm_bandwidth = config.corrector_bandwidth;
    actuactor = tf(1,[1/2/pi/cm_bandwidth 1]);
elseif strcmpi(config.actuator, 'current')
    % Current regulator on corrector magnet (second-order model)
    ps_damping = config.ps_damping;
    ps_omegan = 2*pi*config.corrector_bandwidth/sqrt(1-2*ps_damping^2+sqrt(2-4*ps_damping^2+4*ps_damping^4));
    actuactor = tf(ps_omegan^2,[1 2*ps_damping*ps_omegan ps_omegan^2]);
else
    actuactor = 1;
end

% Vacuum chamber
% First-order model for a thin circular chamber
% More details in Podobedov, B. et al. "Eddy Current Shielding by
% Electrically Thick Vacuum Chamber", PAC'2009

% 0.14e7 % Stainless steel 301 
% 1/2.48e-7 %Steel alloy 4340
% 1/7.2e-7 % Stainless steel 304
% 1/7.4e-7 % Stainless steel 316 
% 1/6e-7 %Stainless steel 405 440A
% 1/8.3e-7 %Stainless steel 17-7PH
% 5.85e7 % Copper 
% 3.45e7 % Aluminum 

mu0 = 4*pi*1e-7;
vac_radius = 11.7e-3;
vac_thickness = 0.1e-3;
vac_conductivity = 5.85e7;
tau = 0.5*mu0*vac_conductivity*vac_radius*vac_thickness;

vacuum_chamber = tf(1,[1/2/pi/config.vac_cutoff 1]);

% Power supply + orbit corrector magnet + vacuum chamber + delays ensemble
corr = c2d(tf(1,1,'ioDelay',latency)*actuactor*vacuum_chamber, Ts, 'zoh');