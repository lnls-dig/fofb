config = struct(...
    'bpm_cic_nstages', 3, ...
    'bpm_cic_ndelays', 1, ...
    'bpm_cic_decimation_factor', 1250, ...
    'actuator', 'current', ...
    'ps_damping', 1, ...
    'corrector_bandwidth', 5e3, ...
    'vac_cutoff', 3.5e3, ...
    'communication_latency', 10e-6, ...
    'ctrl_sampling_rate', 100e3);

[sys_bpm, sys_corr] = fofb_sisomodel(config);

p = struct('actuator_type', 'current', ...
           'disturbance_amplification', 1.7, ...
           'phase_margin', 60, ...
           'communication_latency', config.communication_latency, ...
           'ctrl_sampling_rate', config.ctrl_sampling_rate, ...
           'corrector_bandwidth', config.corrector_bandwidth);
       
r = fofb_bestctrl(p);

controller = r(1,1,1).best_controllers;

result = fofb_freqresp(sys_bpm, sys_corr, controller);