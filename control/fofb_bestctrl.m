function r = fofb_bestctrl(p)

actuator_type = p.actuator_type;
disturbance_amplification = p.disturbance_amplification;
phase_margin = p.phase_margin;
communication_latency = p.communication_latency;
ctrl_sampling_rate = p.ctrl_sampling_rate;
corrector_bandwidth = p.corrector_bandwidth;
vac_cutoff = 3500;
ps_damping = 1;
cic_nstages = 4;
cic_ndelays = 2;
cic_decimation_factor= 1250;

if strcmpi(actuator_type, 'current')
    control_type = 'pid';
elseif strcmpi(actuator_type, 'voltage')
    control_type = 'imc';
end

best_controllers = cell(length(ctrl_sampling_rate), length(communication_latency), length(corrector_bandwidth));
best_results = cell(length(ctrl_sampling_rate), length(communication_latency), length(corrector_bandwidth));
best_crossover_frequencies = zeros(length(ctrl_sampling_rate), length(communication_latency), length(corrector_bandwidth));

crossover_frequency_step = 10;

for i=1:length(ctrl_sampling_rate)
    for j=1:length(communication_latency)
        for k=1:length(corrector_bandwidth)
            config = struct(...
                'bpm_cic_nstages', cic_nstages, ...
                'bpm_cic_ndelays', cic_ndelays, ...
                'bpm_cic_decimation_factor', cic_decimation_factor, ...
                'actuator', actuator_type, ...
                'ps_damping', ps_damping, ...
                'corrector_bandwidth', corrector_bandwidth(k), ...
                'vac_cutoff', vac_cutoff, ...
                'communication_latency', communication_latency(j), ...
                'ctrl_sampling_rate', ctrl_sampling_rate(i));

            [sys_bpm, sys_corr] = fofb_sisomodel(config);
            plant = sys_bpm*sys_corr;

            fprintf('FOFB rate: %0.0f Hz / Comm. delay: %0.0f us / Corrector BW: %0.0f Hz \n', ctrl_sampling_rate(i), 1e6*communication_latency(j), corrector_bandwidth(k));
            crossover_frequency = ctrl_sampling_rate(i)/100;
            while 1
                fprintf('.');
                if strcmpi(control_type, 'pid')
                    pidopt = pidtuneOptions('CrossoverFrequency', 2*pi*crossover_frequency, 'PhaseMargin', phase_margin+5);
                    controller = pidtune(plant, 'PI', pidopt);
                elseif strcmpi(control_type, 'imc')
                    controller = utTuningIMC(plant, 1/(2*pi*crossover_frequency));
                else
                    controller = 1;
                end

                r = fofb_freqresp(sys_bpm, sys_corr, controller);
                
                if strcmpi(r.stability, 'instable') || (r.max_disturbance_amplification > disturbance_amplification) || (r.phase_margin < phase_margin)
                   if ~isempty(best_controllers{i,j,k})
                       fprintf(' Crossover at %0.0f Hz\n\n', crossover_frequency);
                       break;                       
                   else
                       crossover_frequency = crossover_frequency - crossover_frequency_step;
                   end
                else
                    best_controllers{i,j,k} = controller;
                    best_results{i,j,k} = r;
                    best_crossover_frequencies(i,j,k) = crossover_frequency;
                    crossover_frequency = crossover_frequency + crossover_frequency_step;
                end
            end
        end
    end
end

r = struct('best_controllers', best_controllers, ...
           'best_results', best_results, ...
           'best_crossover_frequencies', best_crossover_frequencies);
