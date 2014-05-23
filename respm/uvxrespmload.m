function [matrix, beam_energy, row_names, column_names] = uvxrespmload(filename)

file_id = fopen(filename);

beam_energy = textscan(fgetl(file_id), '%*s%f', 'delimiter', '\t');
beam_energy = beam_energy{1};

column_names = textscan(fgetl(file_id), '%s', 'delimiter', '\t');
column_names = column_names{1};

% Remove empty string of first cell array element
column_names = column_names(2:end);

% Store file position indicator just after first line of numeric data
position = ftell(file_id);

try
    % With the number of columns, build the format string for data
    formatstr = '%*s';
    for i=1:length(column_names)
        formatstr = [formatstr '%f'];
    end

    % Read matrix data
    matrix = textscan(file_id, formatstr, 'delimiter', '\t', 'CollectOutput', 1);
    matrix = matrix{1};

    % Recover file position indicator just after first line
    frewind(file_id);
    fseek(file_id,position,0);

    % Read row names
    row_names = textscan(file_id, '%s%*[^\n]', 'delimiter', '\t', 'CollectOutput', 1);
    row_names = row_names{1};
end

fclose(file_id);
