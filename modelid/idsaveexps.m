function idsaveexps(path, exps, offset)

dateformat = 'yyyy/mm/dd HH:MM:SS';

for i = 1:length(exps)
    start_date = exps(i).start_date;
    stop_date = datestr(datenum(start_date, dateformat) + 50*(exps(i).npts_exp + exps(i).npts_interval)*320e-6/60/60/24  + offset, dateformat);
    start_date = datestr(datenum(start_date, dateformat) - 32/60/60/24  + offset, dateformat);
    
    f = fafind(path, start_date, stop_date);
    data = faload(f, [], [], 100e3);
    
    start_index = find((diff(data.corr_setpoints(:,1)) ~= 0) & (sum(abs(diff(data.corr_setpoints(:,2:end))),2) == 0));
    start_index = start_index(1) + 1 - exps(i).npts_exp - exps(i).npts_interval;
    
    [data_array, corr_iddata, bpm_iddata] = fafcexps(data, start_index, exps(i).npts_exp, exps(i).npts_interval);
    
    dstr = datestr(fatimelvrt2m(data.time(1)), dateformat);
    dstr([findstr(dstr,'/') findstr(dstr,':')]) = '.';
    dstr(findstr(dstr,' ')) = '_';

    save([dstr '.mat'], 'corr_iddata', 'bpm_iddata');
end