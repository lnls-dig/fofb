function fa_data = faloaddata(filenames)
%
% FALOADDATA Loads acquisition data from file.
% fa_data = faloaddata(filename)

if ischar(filenames)
    try
        fileid = fopen(filenames);
        period = fread(fileid,1,'uint32','l');
        header_bpm = strread(fgetl(fileid), '%s', 'delimiter', '\t');
        length_bpm = length(header_bpm);
        header_ps = strread(fgetl(fileid), '%s', 'delimiter', '\t');
        length_ps = length(header_ps);

        data = [];
        while true
            nrows = fread(fileid,1, 'uint32=>double', 'l');
            ncols = fread(fileid,1, 'uint32=>double', 'l');
            if isempty(ncols*nrows)
                break;
            end
            subdata = fread(fileid, ncols*nrows, 'double', 'l');
            data = [data; reshape(subdata, ncols, nrows)'];
        end

        time = typecast(data(:,end), 'uint64');
        data = data(:, 1:end-1);

    catch err
        fclose(fileid);
        rethrow(err);
    end

    fclose(fileid);

    fa_data = struct('time', time, ...
               'bpm_readings', data(:,1:length_bpm), ...
               'bpm_names', {header_bpm(1:end)}, ...
               'ps_readings', data(:,length_bpm+1:length_bpm+length_ps/2), ...
               'ps_names', {header_ps(1:length_ps/2)}, ...
               'ps_setpoints', data(:,1+length_bpm+length_ps/2:length_bpm+length_ps), ...
               'ps_setpoints_names', {header_ps(1+length_ps/2:length_ps)}, ...
               'period', period);

elseif iscell(filenames)
    fa_data = struct('time', [], ...
           'bpm_readings', [], ...
           'bpm_names', [], ...
           'ps_readings', [], ...
           'ps_names', [], ...
           'ps_setpoints', [], ...
           'ps_setpoints_names', [], ...
           'period', []);

    for i=1:length(filenames)
        sub_fa_data = faloaddata(filenames{i});
        
        fa_data.time         = [fa_data.time;         sub_fa_data.time];
        fa_data.bpm_readings = [fa_data.bpm_readings; sub_fa_data.bpm_readings];
        fa_data.ps_readings  = [fa_data.ps_readings;  sub_fa_data.ps_readings];
        fa_data.ps_setpoints = [fa_data.ps_setpoints; sub_fa_data.ps_setpoints];
        
        fa_data.bpm_names = compare(fa_data.bpm_names, sub_fa_data.bpm_names);
        fa_data.ps_names = compare(fa_data.ps_names, sub_fa_data.ps_names);
        fa_data.ps_setpoints_names = compare(fa_data.ps_setpoints_names, sub_fa_data.ps_setpoints_names);
        fa_data.period = compare(fa_data.period, sub_fa_data.period);
    end
end

function new = compare(current, new)

if ~(isempty(current) || isequal(current, new))
    error('Corresponding headers and constants across different files shall be identical.');
end