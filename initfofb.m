cd ..
currentpath = pwd;
cd fofb

pathlist = {...
    fullfile('libdsp'), ...
    fullfile('libunits'), ...
    };

for i = 1:length(pathlist)
    addpath(fullfile(currentpath, pathlist{i}))
end

addpath(genpath(fullfile(currentpath, 'fofb')))
rmpath(genpath(fullfile(currentpath, 'fofb','.git')))
rmpath(genpath(fullfile(currentpath, 'fofb','ofbsim','slprj')))
