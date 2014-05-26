function fadata = faload(filenames)
%
% FALOAD Loads acquisition data from file.
% fadata = faload(filenames)

if ischar(filenames)
    try
        fileid = fopen(filenames);
        period = fread(fileid,1,'uint32','l');
        
        header_bpm = fgetl(fileid);
        if ~isempty(header_bpm)
            header_bpm = textscan(header_bpm, '%s', 'delimiter', '\t');
            header_bpm = header_bpm{1};
        else
            header_bpm = {};
        end        
        length_bpm = length(header_bpm);
        
        header_ps = fgetl(fileid);
        if ~isempty(header_ps)
            header_ps = textscan(header_ps, '%s', 'delimiter', '\t');
            header_ps = header_ps{1};
        else
            header_ps = {};
        end
        length_ps = length(header_ps);
        
        nvars = fread(fileid,1, 'uint32', 'l');

        fileinfo = dir(filenames);
        nrows = (fileinfo.bytes-ftell(fileid))/4/(nvars+2);
        
        data = fread(fileid, nrows*(nvars+2), 'single=>single', 'l');
        data = reshape(data, nvars+2, nrows)';

        time_hi = typecast(data(:,end), 'uint32');
        time_lo = typecast(data(:,end-1), 'uint32');
        time = bitor(bitshift(uint64(time_hi), 32), uint64(time_lo));
                
        data = data(:, 1:end-2);

    catch err
        fclose(fileid);
        rethrow(err);
    end

    fclose(fileid);

    fadata = struct('time', time, ...
               'bpm_readings', data(:,1:length_bpm), ...
               'bpm_names', {header_bpm(1:end)}, ...
               'ps_readings', data(:,length_bpm+1:length_bpm+length_ps/2), ...
               'ps_names', {header_ps(1:length_ps/2)}, ...
               'ps_setpoints', data(:,1+length_bpm+length_ps/2:length_bpm+length_ps), ...
               'ps_setpoints_names', {header_ps(1+length_ps/2:length_ps)}, ...
               'period', period);

elseif iscell(filenames)
    fadata = struct('time', [], ...
           'bpm_readings', [], ...
           'bpm_names', [], ...
           'ps_readings', [], ...
           'ps_names', [], ...
           'ps_setpoints', [], ...
           'ps_setpoints_names', [], ...
           'period', []);

    for i=1:length(filenames)
        subfadata = faload(filenames{i});
        
        fadata.time         = [fadata.time;         subfadata.time];
        fadata.bpm_readings = [fadata.bpm_readings; subfadata.bpm_readings];
        fadata.ps_readings  = [fadata.ps_readings;  subfadata.ps_readings];
        fadata.ps_setpoints = [fadata.ps_setpoints; subfadata.ps_setpoints];
        
        fadata.bpm_names = compare(fadata.bpm_names, subfadata.bpm_names);
        fadata.ps_names = compare(fadata.ps_names, subfadata.ps_names);
        fadata.ps_setpoints_names = compare(fadata.ps_setpoints_names, subfadata.ps_setpoints_names);
        fadata.period = compare(fadata.period, subfadata.period);
    end
end

function new = compare(current, new)

if ~(isempty(current) || isequal(current, new))
    error('Corresponding headers and constants across different files shall be identical.');
end