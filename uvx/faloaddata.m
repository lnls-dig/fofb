function fadata = faloaddata(filenames)
%
% FALOADDATA Loads acquisition data from file.
% fadata = faloaddata(filenames)

if ischar(filenames)
    try
        fileid = fopen(filenames);
        period = fread(fileid,1,'uint32','l');
        header_bpm = textscan(fgetl(fileid), '%s', 'delimiter', '\t');
        header_bpm = header_bpm{1};
        length_bpm = length(header_bpm);
        header_ps = textscan(fgetl(fileid), '%s', 'delimiter', '\t');
        header_ps = header_ps{1};
        length_ps = length(header_ps);

        data = [];
        while true
            nrows = fread(fileid,1, 'uint32=>double', 'l');
            ncols = fread(fileid,1, 'uint32=>double', 'l');
            if isempty(ncols*nrows)
                break;
            end
            subdata = single(fread(fileid, ncols*nrows, 'double', 'l')); % FIXME: remove single() in the future, when data is already encoded in single precision
            %subdata = fread(fileid, ncols*nrows, 'double', 'l');
            data = [data; reshape(subdata, ncols, nrows)'];
        end

        time = typecast(data(:,end), 'uint64');
        data = data(:, 1:end-1);

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
        sub_fadata = faloaddata(filenames{i});
        
        fadata.time         = [fadata.time;         sub_fadata.time];
        fadata.bpm_readings = [fadata.bpm_readings; sub_fadata.bpm_readings];
        fadata.ps_readings  = [fadata.ps_readings;  sub_fadata.ps_readings];
        fadata.ps_setpoints = [fadata.ps_setpoints; sub_fadata.ps_setpoints];
        
        fadata.bpm_names = compare(fadata.bpm_names, sub_fadata.bpm_names);
        fadata.ps_names = compare(fadata.ps_names, sub_fadata.ps_names);
        fadata.ps_setpoints_names = compare(fadata.ps_setpoints_names, sub_fadata.ps_setpoints_names);
        fadata.period = compare(fadata.period, sub_fadata.period);
    end
end

function new = compare(current, new)

if ~(isempty(current) || isequal(current, new))
    error('Corresponding headers and constants across different files shall be identical.');
end