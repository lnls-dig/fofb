%% Parameters
plot_options.naive_bode = false;
plot_options.fit_pct = false;
plot_options.sysid_plot_pause = 0;

nseg_discard_pre = 50;
fir_remove_switching = dfilt.dffir(ones(1,4)/4);
iir_lowpass = design(fdesign.lowpass(11e3, 12e3, 1, 80, Fs), 'butter', 'MatchExactly', 'stopband');

%% Load data
tic
fprintf('Loading data...');
[y, u, Ts, ychname, uchname, seglen] = load_fofbsysid_data_may2023('~/repos/fofb/sysid/mat');
fprintf('DONE. Elapsed time: %f s.\n', toc);

%% Data preprocessing
tic
fprintf('Preprocessing data...\n');
nchan = size(y,2);

%[y, nseg_discard_post] = prbs_remdrift(y, seglen);
nseg_discard_post = 0;

ndiscard = [nseg_discard_pre nseg_discard_post]*seglen;
u = sigcond(u, seglen, iir_lowpass, fir_remove_switching, ndiscard);
y = sigcond(y, seglen, iir_lowpass, fir_remove_switching, ndiscard);
fprintf('DONE. Elapsed time: %f s.\n', toc);

%% System Identification
sysid_result = {};
sys = {};
fitpct = zeros(nchan,1);

arx_orders = repmat([6 6 2], nchan, 1);
arx_orders(87, :) =  [14 14 2];
arx_orders(134, :) =  [6 6 1];

for i=1:nchan
    tic
    fprintf('Estimating %s -> BPM model (idx = %4d)... ', uchname{i}, i);
    sysid_data = iddata(detrend(y(:,i),0), detrend(u(:,i),0), Ts, 'OutputName', ychname{i}, 'InputName', uchname{i});
    sysid_result{i} = arx(sysid_data, [arx_orders(i, 1:2) delayest(sysid_data)], arxOptions('Focus', 'Simulation', 'EnforceStability', true));
    fitpct(i) = sysid_result{i}.Report.Fit.FitPercent;
    fprintf('Fit = %6.2f%%, Elapsed time: %f s.\n', fitpct(i), toc);    
    sys{i} = ss(sysid_result{i});    
end

%% Normalize DC gain to 1
for i=nchan:-1:1
    sys_norm{i} = sys{i}/dcgain(sys{i});
end

%% Plot results

%
if plot_options.sysid_plot_pause > 0
    for i=1:nchan        
        %figure
        compare(sysid_data, sysid_result{i})
        pause(plot_options.sysid_plot_pause)
    end
end

% Naive Bode plot
if plot_options.naive_bode
    Y_avg = fft(y_avg);
    U_avg = fft(u_avg);
    
    Y_avg = Y_avg(1:seglen/2+1, :);
    U_avg = U_avg(1:seglen/2+1, :);
    
    f = (0:seglen/2)'*Fs/(seglen+1);
    
    amplitude_norm = Y_avg(2,:)./U_avg(2,:);
    
    figure
    ax(1) = subplot(211);
    plot(f, 20*log10(abs(Y_avg./U_avg./amplitude_norm)))
    ylabel('Magnitude [dB]')
    grid on
    ax(2) = subplot(212);
    plot(f, 180/pi*unwrap(angle(Y_avg./U_avg)));
    ylabel('Phase [deg]')
    grid on
    xlabel('Frequency [Hz]')
    linkaxes(ax,'x')
end

if plot_options.fit_pct
end