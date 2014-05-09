function haxis_out = plot_spectrum(signals, Fs, graph_type, amplitude_range, freq_range, window, color, log, signal_names, unit, dataset_name, haxis_in)
% haxis_out = fofb_plot_spectrum(signals, Fs, graph_type, amplitude_range, freq_range, window, color, log, signal_names, unit, dataset_name, haxis_in)

[npts, nsignals] = size(signals);

if nargin < 2 || isempty(Fs)
    Fs = 1;
end

if nargin < 3 || isempty(graph_type)
    graph_type = 'dft';
end

if (nargin < 5) || isempty(freq_range)
    freq_range = [0 Fs/2];
end

if (nargin < 6) || isempty(window)
    window = [];
else
    PG_normalization = sum(window)/npts;
    window = window/PG_normalization;
end

if (nargin < 7) || isempty(color)
    color = 'b';
end

if isscalar(color) || (size(color,1) ~= nsignals)
    color = repmat(color, nsignals, 1);
end

if (nargin < 8) || isempty(log)
    log = 1;
end

if (nargin < 9)
    signal_names = [];
end

if (nargin < 10) || isempty(unit)
    unit = 'a.u.';
end

if (nargin < 11) || isempty(unit)
    dataset_name = [];
end

% ============
% Calculations
% ============


if (nargin < 1) || isempty(signals)
    spectra = [];
    freq = [];
elseif strcmpi(graph_type, 'dft')
    [spectra, freq] = fourierseries(signals, Fs, window);
    spectrum_ylabel = sprintf('Amplitude (%s)', unit);
elseif strcmpi(graph_type, 'psd')
    [spectra, freq] = psdrms(signals, Fs, 0, Inf, window, [], [], graph_type);
    spectrum_ylabel = sprintf('PSD (%s/\\surdHz)', unit);
elseif strcmpi(graph_type, 'rms')
    [spectra, freq] = psdrms(signals, Fs, 0, Inf, window, [], [], graph_type);
    spectrum_ylabel = sprintf('Integrated RMS (%s)', unit);
else
    error('Invalid value for ''graph_type'' parameter.');
end

if (nargin < 4) || isempty(amplitude_range)
    amplitude_range = [min(min(spectra(2:end,:))) max(max(spectra))];
end

% ===========
% Plot graphs
% ===========

if (nargin < 12) || isempty(haxis_in)
    haxis = zeros(nsignals, 1);
else
    haxis = haxis_in;
end

for i = 1:nsignals
    if (nargin < 12) || isempty(haxis_in)
        figure;
        haxis(i) = gca;
    else
        axes(haxis(i));
    end
    if ~isempty(spectra)
        hnew = plot(freq, spectra(:,i), 'Color', color(i, :));
    end    
    hold on;
    grid on;
    if log
        set(haxis(i), 'YScale', 'log');
    else
        set(haxis(i), 'YScale', 'linear');
    end
    axis([freq_range amplitude_range])
    ylabel(spectrum_ylabel,'FontSize',14,'FontWeight','bold');
    xlabel('Frequency (Hz)','FontSize',14,'FontWeight','bold');    
    set(haxis(i), 'FontSize', 14);
    if ~isempty(signal_names) && ~isempty(spectra)
        [~,~,outh,outm] = legend;
        n = length(outm);
        outm{n+1} = '';
        if ~isempty(signal_names)
            outm{n+1} = signal_names{i};
        end
        if ~isempty(dataset_name)
            outm{n+1} = [outm{n+1} ' - ' dataset_name];
        end
        
        legend([outh;hnew],outm, 'FontSize', 10);
    end
end

if nargout > 0
    haxis_out = haxis;
end