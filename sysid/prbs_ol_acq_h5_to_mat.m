% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc; clear;

prbs_ol_acq_fpath = 'current/';
ps_names_fpath = 'ps_names.txt';
ps_names_text = fileread(ps_names_fpath);
ps_names_list = regexp(ps_names_text,'[^\n]*[^\n]','match')';
ncorr = size(ps_names_list,1);

excluded_corr = [1, 80 ,81 ,160];

for i=1:ncorr
    fprintf('Corrector %d\n',i);
    tic;
    ps_name = ps_names_list{i};
    fpath = [prbs_ol_acq_fpath ps_name '.h5'];
    if ~ismember(i, excluded_corr)
        data.prbs_lfsr_len = h5read(fpath,'/data/prbs_lfsr_len');
        data.prbs_step_duration = h5read(fpath,'/data/prbs_step_duration');
        data.prbs_mov_avg_taps = h5read(fpath,'/params/prbs_mov_avg_taps');
        data.sampling_frequency = h5read(fpath,'/data/sampling_frequency');

        acq_nrpoints_before = h5read(fpath,'/params/acq_nrpoints_before');

        % PRBS input data
        prbs_data = double(h5read(fpath,'/data/prbs_data/0'));
        prbs_data = prbs_data(acq_nrpoints_before+1:end);
        corr_type_label = fpath(6:7);
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
        orbx = h5read(fpath,'/data/orbx');
        orby = h5read(fpath,'/data/orby');
        data.orbx = orbx(:,acq_nrpoints_before+1:end);
        data.orby = orby(:,acq_nrpoints_before+1:end);

        save(ps_name,'data');
    end
    fprintf('Elapsed time: %f s\n',toc);
end
