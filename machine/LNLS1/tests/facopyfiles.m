function facopyfiles(path, filenames)

for i=1:length(filenames)
    [~, name, ext] = fileparts(filenames{i});
    copyfile(filenames{i}, 'temp.dat');
    movefile('temp.dat', fullfile(path, [name ext]));
end