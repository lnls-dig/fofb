function r = faverify(filenames)

tlast = [];
j = 1;
for i=1:length(filenames)
    [pathstr, filename, ext] = fileparts(filenames{i}) ;
    
    if strcmpi(ext, '.dat') && ~strcmpi(filename, 'temp')
        fa_data = faloaddata(filenames{i});
        t = [tlast; fa_data.time];
        dt = diff(t);
        tlast = fa_data.time(end);

        failed = dt > fa_data.period*1e3*1.5;
        if any(failed)
            r(j).filename = filenames{i};
            r(j).failed_transition = (failed(1) == 1);
            r(j).failed = find(failed(2:end))+1;
            j=j+1;
        end
    end
end
