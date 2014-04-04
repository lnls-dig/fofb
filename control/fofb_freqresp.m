function r = fofb_freqresp(sys_bpm, sys_corr, controller)

w = warning('off', 'Control:ltiobject:UseSSforInternalDelay');
UseSSforInternalDelay_state = w.state;
w = warning('off', 'Control:analysis:MarginUnstable');
MarginUnstable_state = w.state;

plant = sys_bpm*sys_corr;
open_loop = controller*plant;
closed_loop_beam = feedback(controller*sys_corr, sys_bpm);
sensitivity = feedback(1, controller*plant);

Ts = get(sys_bpm, 'Ts');

% Define frequency range to be evaluated and frequency range for plotting
frequency = logspace(log10(min(0.01, 1/Ts/2/2/1000)), log10(1/Ts/2), 1000)';

% Calculate open-loop frequency response (magnitude and phase)
[mag_open_loop, phase_open_loop] = bode(open_loop, 2*pi*frequency);
mag_open_loop = squeeze(mag_open_loop);
phase_open_loop = squeeze(phase_open_loop);
[Gm,Pm,Wcg,Wcp] = margin(open_loop);

% Calculate closed-loop frequency response (magnitude)
mag_closed_loop = squeeze(bode(closed_loop_beam, 2*pi*frequency));

% Calculate disturbance rejection transfer function (magnitude)
mag_disturbance_rejection = squeeze(bode(sensitivity, 2*pi*frequency));

% Find maximum disturbance amplification
[max_disturbance_amplification, index_max_disturbance_amplification] = max(mag_disturbance_rejection);
index_disturbance_amplification = find(mag_disturbance_rejection > 1);

% Find disturbace attenuation-amplification boundaries (frequencis for each the disturbance rejection, in dB, changes its signal)
crossover_points = index_disturbance_amplification((index_disturbance_amplification(2:end)-index_disturbance_amplification(1:end-1))~=1);
if ~isempty(index_disturbance_amplification)
    if index_disturbance_amplification(1) > 1
        crossover_points = [index_disturbance_amplification(1); crossover_points];
    end
    if index_disturbance_amplification(end) < length(mag_disturbance_rejection)
        crossover_points = [crossover_points; index_disturbance_amplification(end)];
    end
end

if isempty(crossover_points)
    crossover_points = 1;
end

% Check closed-loop stability
if isstable(closed_loop_beam)
    closed_loop_stability = 'stable';
else
    closed_loop_stability = 'instable';
end

r = struct('frequency', frequency, ...
           'mag_disturbance_rejection', mag_disturbance_rejection, ...
           'phase_margin', Pm, ...
           'stability', closed_loop_stability, ...
           'crossover_frequency', frequency(crossover_points(1)), ...
           'max_disturbance_amplification', max_disturbance_amplification);            

if nargout == 0
    frequency_plot_range = [min(0.1, 1/Ts/2/2/100) 1/Ts/2];
    
    N = 2;
    if length(crossover_points) > N
        crossover_points = crossover_points(1:N);
    end
    
    % Plot results
    screensize = get(0, 'ScreenSize');
    h = figure('Position', [0 0 screensize(3) screensize(4)]);
    set(h, 'Name', 'Fast Orbit Feedback Simulator - System frequency response');

    % Plot disturbance rejection
    ax = subplot(311);
    semilogx(frequency, 20*log10(mag_disturbance_rejection), 'r');
    grid('on');
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    legend('Disturbance suppression');
    hold on;
    text(frequency(index_max_disturbance_amplification), 20*log10(max_disturbance_amplification), ['Maximum amplification: ' num2str(max_disturbance_amplification) ' (' num2str(20*log10(max_disturbance_amplification)) ' dB)  '], 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right', 'Color', 'r', 'FontSize', 8, 'Clipping', 'on');
    text(frequency(index_max_disturbance_amplification), 0, [num2str(frequency(index_max_disturbance_amplification)) ' Hz  '], 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'r', 'FontSize', 8, 'Clipping', 'on');
    stem(frequency(index_max_disturbance_amplification), 20*log10(max_disturbance_amplification), 'ro', 'filled');
    for i=1:length(crossover_points)
        freq = frequency(crossover_points(i));
        text(freq, 20*log10(mag_disturbance_rejection(crossover_points(i))), [num2str(frequency(crossover_points(i))) ' Hz  '], 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'r', 'FontSize', 8, 'Clipping', 'on');
    end
    plot(frequency(crossover_points), 20*log10(mag_disturbance_rejection(crossover_points)), 'r*');    
    set(ax, 'FontSize', 8, 'XLim', frequency_plot_range, 'YLim', [-80 20*ceil(log10(max_disturbance_amplification))]);

    % Plot open- and closed-loop gains
    ax = subplot(312);
    semilogx(frequency, 20*log10([mag_open_loop mag_closed_loop]));
    grid('on');
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    h_legend = legend('Open loop', ['Closed loop (' closed_loop_stability ')']);
    if ~isstable(closed_loop_beam)
        set(h_legend, 'Color', [1 0.9 0.9]);
        set(h_legend, 'EdgeColor', [1 0 0]);
    end
    set(ax, 'FontSize', 8, 'XLim', frequency_plot_range, 'YLim', [-150 150]);
    hold on;
    text(Wcg/2/pi, -20*log10(Gm), ['Gain margin: ' num2str(Gm) ' (' num2str(20*log10(Gm)) ' dB)  '], 'VerticalAlignment', 'top', 'HorizontalAlignment','right', 'Color', 'b', 'FontSize', 8, 'Clipping', 'on');
    stem(Wcg/2/pi, 20*log10(interp1(frequency,mag_open_loop,Wcg/2/pi)), 'bo', 'filled');

    % Plot open-loop phase
    ax = subplot(313);
    semilogx(frequency, phase_open_loop);
    grid('on');
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)');
    legend('Open loop');
    set(ax, 'FontSize', 8, 'XLim', frequency_plot_range, 'YLim', [-1000 1000]);
    hold on;
    text(Wcp/2/pi, interp1(frequency,phase_open_loop,Wcp/2/pi), ['  Phase margin: ' num2str(Pm) '°'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left', 'Color', [0 0 1], 'FontSize', 8, 'Clipping', 'on');
    plot(Wcp/2/pi, interp1(frequency,phase_open_loop,Wcp/2/pi), 'o', 'MarkerFaceColor', 'b');
end

warning(UseSSforInternalDelay_state, 'Control:ltiobject:UseSSforInternalDelay');
warning(MarginUnstable_state, 'Control:analysis:MarginUnstable');