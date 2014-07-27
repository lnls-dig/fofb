function idsaveexps(path, start_date_array, npts_exp, npts_interval, offset)

dateformat = 'yyyy/mm/dd HH:MM:SS';

for i = 1:length(start_date_array)
    start_date = start_date_array{i};
    stop_date = datestr(datenum(start_date, dateformat) + 3/60/24  + offset, dateformat);
    start_date = datestr(datenum(start_date, dateformat) - 32/60/60/24  + offset, dateformat);
    
    f = fafind(path, start_date, stop_date);
    data = faload(f, [], [], 100e3);
    
    start_index = find((diff(data.corr_setpoints(:,1)) ~= 0) & (sum(abs(diff(data.corr_setpoints(:,2:end))),2) == 0));
    start_index = start_index(1) + 1 - npts_exp - npts_interval;
    
    [data_array, corr_iddata, bpm_iddata] = fafcexps(data, start_index, npts_exp, npts_interval);
    
    dstr = start_date_array{i};
    dstr([findstr(dstr,'/') findstr(dstr,':')]) = '.';
    dstr(findstr(dstr,' ')) = '_';

    save([dstr '.mat'], 'corr_iddata', 'bpm_iddata');
end