% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>
% Modified by: Guilherme Ricioli  <guilherme.ricioli@lnls.br>

clc; clear;

%acq_fpath = 'rf-phase-fofb-off.h5';
acq_fpath = '/ibira/lnls/labs/swc/MachineStudies/22-07-2024/rf-phase-fofb-off.h5';
% NOTE: This has to be a dispersive BPM
bpm = 19;

orbx = h5read(acq_fpath,'/data/orbx');
sampling_frequency = h5read(acq_fpath,'/data/sampling_frequency');
Ts = 1/sampling_frequency;

% Model for fitting
fun = @(x,xdata) ...
      x(1).*exp(-x(2).*xdata).*sin(x(3).*xdata) +...
      x(4).*exp(-x(2).*xdata).*cos(x(3).*xdata);
x0 = [100 1300 2*pi*2000 100];

% Filter for switching noise removal
mov_avg_filt = dfilt.dffir(ones(1,4)/4);

% Step response windows visually chosen
step_resp_ivals = {22341:22473;
                   124193:124350;
                   172580:172721;
                   177210:177340;
                   181840:181960;
                  };
n_trials = size(step_resp_ivals,1);

orbx = filter(mov_avg_filt,double(orbx(bpm,:)));

alpha_s_avg = 0;
big_omega_avg = 0;
for i=1:n_trials
    fprintf('Trial %d: ', i);

    y = orbx(step_resp_ivals{i});
    y = detrend(y,0);
    t = linspace(0,length(y)*Ts,length(y));

    % Performs fit
    x = lsqcurvefit(fun,x0,t,y);

    fprintf('Omega: %.2f Hz, alpha_s: %.2f s^-1\n', x(3)/(2*pi), x(2));

    figure;
    title(string(i));
    hold on;
    scatter(t,y);
    plot(t,fun(x,t));
    xlabel('Time [s]');
    ylabel('Position displacement [um]')
    legend({'Data','Fit'});
    grid on;

    alpha_s_avg = alpha_s_avg + x(2)/n_trials;
    big_omega_avg = big_omega_avg + x(3)/n_trials;
end

fprintf('\n');
fprintf('Omega (average): %.2f Hz\n', big_omega_avg/(2*pi));
fprintf('alpha_s (average): %.2f s^-1\n', alpha_s_avg);

save('rf_fitted_params','big_omega_avg','alpha_s_avg');
