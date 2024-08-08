% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc; clear;

prbs_ol_acq_fpath = 'current/';
ps_names_fpath = 'ps_names.txt';
ps_name_bpm_idx_max_response_fpath = 'relation_corrector_to_bpm_idx_max_response.txt';

ps_names_table = readtable(ps_names_fpath, Delimiter=',', ...
                           ReadVariableNames=false, TextType='string');
ps_name_bpm_idx_max_response_table = ...
  readtable(ps_name_bpm_idx_max_response_fpath, Delimiter=',', ...
            ReadVariableNames=false, ReadRowNames=true);
ncorr = size(ps_names_table, 1);

excluded_corr = [1, 80, 81, 160];

for i=1:ncorr
    tic;
    ps_name = string(ps_names_table{i, 1});
    fprintf('Corrector %d: %s\n', i, ps_name);
    fpath = strcat(prbs_ol_acq_fpath, ps_name, '.h5');
    if ~ismember(i, excluded_corr)
        data.prbs_lfsr_len = h5read(fpath, '/data/prbs_lfsr_len');
        data.prbs_step_duration = h5read(fpath, '/data/prbs_step_duration');
        data.prbs_mov_avg_taps = h5read(fpath, '/params/prbs_mov_avg_taps');
        data.sampling_frequency = h5read(fpath, '/data/sampling_frequency');

        acq_nrpoints_before = h5read(fpath, '/params/acq_nrpoints_before');

        % PRBS input data
        prbs_data = double(h5read(fpath, '/data/prbs_data/0'));
        prbs_data = prbs_data(acq_nrpoints_before+1:end);
        corr_type_label = ps_name{1}(6:7);
        for j=1:length(prbs_data)
            if prbs_data(j)==0
                prbs_data(j)=-1;
            end
        end
        prbs_data = prbs_data.*20e-3;
        if strcmp(corr_type_label,'C3')
            prbs_data = prbs_data./2;
        end
        data.prbs_data = prbs_data;

        % Orbit output data
        orbx = h5read(fpath, '/data/orbx');
        orby = h5read(fpath, '/data/orby');
        data.orb = [orbx(:, acq_nrpoints_before+1:end)
                    orby(:, acq_nrpoints_before+1:end)];

        data.bpm_idx_max_response = ...
          ps_name_bpm_idx_max_response_table{[ps_name], 1} + 1;

        save(ps_name, 'data');
    end
    fprintf('Elapsed time: %f s\n', toc);
end
