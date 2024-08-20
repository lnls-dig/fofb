% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>
% Modified by: Guilherme Ricioli  <guilherme.ricioli@lnls.br>

clc; clear;

%% Loads MATLAB objects
fofb_mdl_fpath = 'fofb_mdl.mat';
exp_sigma_sim_res_fpath.open_loop = 'exp_sigma_sim-Open Loop.mat';
exp_sigma_sim_res_fpath.sensitivity = 'exp_sigma_sim-Sensitivity.mat';
exp_sigma_res_fpath.open_loop = 'exp_sigma_f-open-loop.mat';
exp_sigma_res_fpath.sensitivity = 'exp_sigma_f-sensitivity.mat';

fofb_mdl = load(fofb_mdl_fpath).fofb_mdl;
exp_sigma_sim_res.open_loop = load(exp_sigma_sim_res_fpath.open_loop);
exp_sigma_sim_res.sensitivity = load(exp_sigma_sim_res_fpath.sensitivity);
exp_sigma_res.open_loop = load(exp_sigma_res_fpath.open_loop);
exp_sigma_res.sensitivity = load(exp_sigma_res_fpath.sensitivity);

n_modes = min(size(fofb_mdl.G));

%% Plots open loop results
figure;

n_plots = 0;

% MATLAB
h = sigmaplot(fofb_mdl.G);
plotoptions = sigmaoptions;
plotoptions.FreqUnits = 'Hz';
plotoptions.Grid = 'on';
setoptions(h,plotoptions);
hold on;

n_plots = n_plots + 1;
leg{n_plots} = 'Sigma (MATLAB)';

% Simulated
for k=1:n_modes
    % We removed DC from PRBS and orbit signals, so start from frequency index 2
    semilogx(exp_sigma_sim_res.open_loop.freqs(2:end),...
             20*log10(squeeze(exp_sigma_sim_res.open_loop.exp_sigma(k,k,2:end))),...
             'Color','#D95319');
    hold on;

    n_plots = n_plots + 1;
    if k==1
      leg{n_plots} = 'Simulated';
    else
      leg{n_plots} = '';
    end
end

% Experimental
% TODO: change sysid/python_lib/fofb_experimental_sigma.ipynb in
%       exp-sigma-notebook branch so that open-loop file doesn't use 'mat_d'
%       struct
for k=1:size(exp_sigma_res.open_loop.mat_d.exp_sigma_f,2)
    semilogx(exp_sigma_res.open_loop.mat_d.freqs,...
             20*log10(exp_sigma_res.open_loop.mat_d.exp_sigma_f(:,k)),...
             'Color','#EDB120');
    hold on;

    n_plots = n_plots + 1;
    if k==1
      leg{n_plots} = 'Experimental';
    else
      leg{n_plots} = '';
    end
end

legend(leg)

% set(gcf,'Units','Inches');
% figpos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [figpos(3), figpos(4)])
% print(gcf,'open_loop_sigma.pdf','-dpdf','-fillpage')

%% Plots sensitivity results
figure;

n_plots = 0;

% MATLAB
h = sigmaplot(fofb_mdl.P('yd','d'));
setoptions(h,plotoptions);
hold on;

n_plots = n_plots + 1;
leg{n_plots} = 'Sigma (MATLAB)';

% Simulated
for k=1:n_modes
    % We removed DC from PRBS and orbit signals, so start from frequency index 2
    semilogx(exp_sigma_sim_res.sensitivity.freqs(2:end),...
             20*log10(squeeze(exp_sigma_sim_res.sensitivity.exp_sigma(k,k,2:end))),...
             'Color','#D95319');

    n_plots = n_plots + 1;
    if k==1
      leg{n_plots} = 'Simulated';
    else
      leg{n_plots} = '';
    end
end

% Experimental
for k=1:size(exp_sigma_res.sensitivity.exp_sigma_f,2)
    semilogx(exp_sigma_res.sensitivity.freqs,...
             20*log10(exp_sigma_res.sensitivity.exp_sigma_f(:,k)),...
             'Color','#EDB120');

    n_plots = n_plots + 1;
    if k==1
      leg{n_plots} = 'Experimental';
    else
      leg{n_plots} = '';
    end
end

legend(leg)

% set(gcf,'Units','Inches');
% figpos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [figpos(3), figpos(4)])
% print(gcf,'sensitivity_sigma.pdf','-dpdf','-fillpage')
