function fadata = faload(filenames, selected_bpm, selected_corr)
%
% FALOAD Loads fast acquisition (FA) data from file.

if nargin < 2
    selected_bpm = [];
end

if nargin < 3
    selected_corr = [];
end

if ischar(filenames)
    try
        fileid = fopen(filenames);
        period = fread(fileid, 1, 'uint32', 'l');

        header_bpm = readheader(fileid);
        nbpm = length(header_bpm);
        
        header_corr = readheader(fileid);
        ncorr = length(header_corr)/2;
        
        if nargin < 2 || isempty(selected_bpm)
            selected_bpm = 1:nbpm;
        elseif selected_bpm == 0
            selected_bpm = [];
        end

        if nargin < 3 || isempty(selected_corr)
            selected_corr = 1:ncorr;
        elseif selected_corr == 0
            selected_corr = [];
        end
        
        header_bpm = header_bpm(selected_bpm);
        header_corr = {header_corr{selected_corr} header_corr{ncorr+selected_corr}}';
        
        nselected_bpm = length(header_bpm);
        nselected_corr = length(header_corr);
        
        data = [];
        while true
            nrows = fread(fileid, 1, 'uint32=>uint32', 'l');
            ncols = fread(fileid, 1, 'uint32=>uint32', 'l');
            if isempty(ncols) || isempty(nrows)
                break;
            end
            subdata = fread(fileid, ncols*nrows, 'single=>single', 'l');
            subdata = reshape(subdata, ncols, nrows)';
            data = [data; subdata(:, [selected_bpm (nbpm+selected_corr) (nbpm+ncorr+selected_corr)]) subdata(:, end-1:end)];
        end

        time_hi = typecast(data(:,end), 'uint32');
        time_lo = typecast(data(:,end-1), 'uint32');
        time = bitor(bitshift(uint64(time_hi), 32), uint64(time_lo));

        data = data(:,1:end-2);

    catch err
        fclose(fileid);
        rethrow(err);
    end

    fclose(fileid);

    fadata = struct('time', time, ...
               'bpm_readings', data(:,1:nselected_bpm), ...
               'bpm_names', {header_bpm(1:end)}, ...
               'corr_readings', data(:,nselected_bpm+1:nselected_bpm+nselected_corr/2), ...
               'corr_names', {header_corr(1:nselected_corr/2)}, ...
               'corr_setpoints', data(:,1+nselected_bpm+nselected_corr/2:nselected_bpm+nselected_corr), ...
               'corr_setpoints_names', {header_corr(1+nselected_corr/2:nselected_corr)}, ...
               'period', period);

elseif iscell(filenames)
    fadata = struct('time', [], ...
           'bpm_readings', [], ...
           'bpm_names', [], ...
           'corr_readings', [], ...
           'corr_names', [], ...
           'corr_setpoints', [], ...
           'corr_setpoints_names', [], ...
           'period', []);

    for i=1:length(filenames)
        sub_fadata = faload(filenames{i}, selected_bpm, selected_corr);

        fadata.time         = [fadata.time;         sub_fadata.time];
        fadata.bpm_readings = [fadata.bpm_readings; sub_fadata.bpm_readings];
        fadata.corr_readings  = [fadata.corr_readings;  sub_fadata.corr_readings];
        fadata.corr_setpoints = [fadata.corr_setpoints; sub_fadata.corr_setpoints];

        fadata.bpm_names = compare(fadata.bpm_names, sub_fadata.bpm_names);
        fadata.corr_names = compare(fadata.corr_names, sub_fadata.corr_names);
        fadata.corr_setpoints_names = compare(fadata.corr_setpoints_names, sub_fadata.corr_setpoints_names);
        fadata.period = compare(fadata.period, sub_fadata.period);
    end
end

function new = compare(current, new)

if ~(isempty(current) || isequal(current, new))
    error('Corresponding headers and constants across different files shall be identical.');
end

function header = readheader(fileid)

header = fgetl(fileid);
if ~isempty(header)
    header = textscan(header, '%s', 'delimiter', '\t');
    header = header{1};
else
    header = {};
end
