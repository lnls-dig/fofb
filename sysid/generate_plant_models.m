% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
ps_names_fpath = "ps_names.txt";
text = fileread(ps_names_fpath);
expr = '[^\n]*[^\n]';
matches = regexp(text,expr,'match')';
ncorr = size(matches,1);

sys=cell(ncorr, 1);
fit = zeros(160,1);

for i=1:ncorr
    fprintf("Corrector %d\n", i);
    tic;
    file=matches{i};
    file = [file '.mat'];
    if (i==1)||(i==80)||(i==81)||(i==160) %excluded correctors 1 80 81 160
        sys{i}=NaN;
    elseif i<= floor(ncorr/2)
        [sys_it,fit_it] = plant_arx_fit(file, [6 6 2], 'X');
        sys{i}=sys_it;
        fit(i)=fit_it;
    elseif i>floor(ncorr/2)
        [sys_it,fit_it] = plant_arx_fit(file, [6 6 2], 'Y');
        sys{i}=sys_it;
        fit(i)=fit_it;
    end
    fprintf("Elapsed time: %f s\n", toc);
end

plot(fit);
title("Fit x corrector number")
xlabel('Corrector number')
ylabel('Fit')

save('sysid_res','sys');