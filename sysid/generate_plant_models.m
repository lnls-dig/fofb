% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>
% Modified by: Guilherme Ricioli <guilherme.ricioli@lnls.br>

clc; clear;

prbs_ol_acq_fpath = '/current/';
ps_names_fpath = 'ps_names.txt';

ps_names_table = readtable(ps_names_fpath, Delimiter=',', ...
                           ReadVariableNames=false, TextType='string');
ncorr = size(ps_names_table, 1);

excluded_corr = [1, 80, 81, 160];

sys = cell(ncorr, 1);
fit = zeros(ncorr, 1);

figure();
for i=1:ncorr
    tic;
    ps_name = string(ps_names_table{i, 1});
    ps_name_without_colon = replace(ps_name, ":", "-");
    fprintf('Corrector %d: %s\n', i, ps_name);
    if ismember(i, excluded_corr)
        sys{i} = NaN;
    else
        fpath = strcat(prbs_ol_acq_fpath, ps_name_without_colon, '.mat');
        sys{i} = fit_plant_arx(fpath, [6 6 2]);
        fit(i) = sys{i}.Report.Fit.FitPercent;

        %opts = bodeoptions;
        %opts.FreqUnits = 'Hz';
        %bode(sys{i}, opts);
        %hold on;
    end
    fprintf('Elapsed time: %f s\n', toc);
end

save('sysid_res', 'sys');

figure();
plot(fit);
title('Fit Percent');
xlabel('Corrector index');
ylabel('Fit [%]');
