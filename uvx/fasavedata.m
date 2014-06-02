function fasavedata(filename, data, timestamp, header, mode)

if nargin < 5 || isempty(mode)
    mode = 'binary';
end

if ischar(filename)
    try
        fileid = fopen(filename, 'w+');
        fprintf(fileid, '%s\t', header{:});
        fprintf(fileid, '\n');

        if strcmpi(mode, 'binary')
            data = [data timestamp];
            dim = uint32(size(data));
            fwrite(fileid, dim, 'uint32', 0, 'l');
            fwrite(fileid, data, 'double', 0, 'l');
        elseif strcmpi(mode, 'ascii')
            for i=1:size(data,1)
                %fprintf(fileid, '%s\t', datestr(timestamp(i), 'yyyy-mm-dd HH:MM:SS.FFFFFF'));
                fprintf(fileid, '%d\t', timestamp(i));
                fprintf(fileid, '%0.10f\t', data(i, :));
                fprintf(fileid, '\n');
            end
        else
            error('Unknown mode.');
        end
    catch err
        fclose(fileid);
        rethrow(err);
    end

    fclose(fileid);
end
