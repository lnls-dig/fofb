function [y, u, Ts, ychname, uchname, seglen] = load_fofbsysid_data_may2023(base_path)

% Prepare filenames
for i=77:-4:1
    temp_names{i,1} = sprintf('%02dM1_FCH', (i+3)/4);
    temp_names{i+1,1} = sprintf('%02dM2_FCH', (i+3)/4);
    temp_names{i+2,1} = sprintf('%02dC2_FCH', (i+3)/4);
    temp_names{i+3,1} = sprintf('%02dC3_FCH', (i+3)/4);
    temp_names{80+i,1} = sprintf('%02dM1_FCV', (i+3)/4);
    temp_names{80+i+1,1} = sprintf('%02dM2_FCV', (i+3)/4);
    temp_names{80+i+2,1} = sprintf('%02dC2_FCV', (i+3)/4);
    temp_names{80+i+3,1} = sprintf('%02dC3_FCV', (i+3)/4);
end

uchname = temp_names;
uchname(1:79) = temp_names(2:80);
uchname(80) = temp_names(1);
uchname(81:159) = temp_names(82:160);
uchname(160) = temp_names(81);

ychname = uchname;

for i=1:length(ychname)
    ychname{i} = ['BPM - ' ychname{i}];
end

% Load data
failed_idx = [];
for i=1:length(uchname)
    try
        data = load(fullfile(base_path, sprintf('%s.mat', uchname{i})));
        data = data.data;
        if strcmpi(uchname{i}(end), 'H')
            y(:,i) = double(data.orbx);
        elseif strcmpi(uchname{i}(end), 'V')
            y(:,i) = double(data.orby);
        end
        u(:,i) = double(data.currdata);
    catch
        failed_idx(end+1) = i;
    end
end

for i=failed_idx
    u(:,i) = 0;
    y(:,i) = 0;
end

lfsr_len = double(data.prbs_lfsr_len(1));
step_duration = double(data.prbs_step_duration(1));

seglen = (2^lfsr_len-1)*step_duration;

Ts = 1/data.sampling_frequency;
