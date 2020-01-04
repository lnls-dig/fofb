clear
%close all

% ==========
% Parameters
% ==========
data_filename = 'vibration_data.mat';
data_unit = 'velocity';

flow = 0;
fhigh = Inf;
wdwfcn = @hamming;
wdwsize = 1/10;  % number of points of each Welch's segment, expressed as a fraction of the number of points of the data vector
overlap = 1/2;   % number of overlapped points between adjacent segments, expressed as fraction of the number of points of the Welch's segment size

armodel_order = 80;
psdestfcnset = {@pyulear, @pburg, @pcov, @pmcov};
arestfcnset = {@aryule, @arburg, @arcov, @armcov};
arlegendset = {'Yule-Walker', 'Burg', 'Covariance', 'Modified covariance'};

% ============
% Process data
% ============

% Load t, vib1 and vib2 variables
load(data_filename);

data = [vib1 vib2];

for i=1:size(data,2)
    x = data(:,i);
    x = x-mean(x);
    x = x/std(x);

    fs = 1/(t(2)-t(1));

    npts = length(x);
    nptswdw = floor(npts*wdwsize);          
    noverlap = floor(nptswdw*overlap);  

    window = wdwfcn(nptswdw);               % window vector

    % Compute PSD
    [Pxx,f] = pwelch(x, window, noverlap, nptswdw, fs, 'onesided');

    % Crop frequency vector at interest bandwidth
    index = find((f >= flow-eps) & (f <= fhigh+eps));
    fsel = f(index);
    df = fsel(2) - fsel(1);
    Pxxsel = Pxx(index);

%     % Compute integrated RMS
%     intrms = sqrt(cumtrapz(Pxxsel)*df);
% 
%     % Compute reverse integrated RMS
%     reversedPxxsel = Pxxsel(end:-1:1);
%     reveresedintrms = sqrt(cumtrapz(reversedPxxsel)*df);
%     reveresedintrms = reveresedintrms(end:-1:1);

    % Estimate AR models and estimate PSD from AR models
    whitenoise = randn(npts,1);
    for j=1:length(psdestfcnset)
        psdestfcn = psdestfcnset{j};
        arestfcn = arestfcnset{j};

        [Pxx_, f_] = psdestfcn(x, armodel_order, nptswdw, fs, 'onesided');
        [arcoeff_, xxxx] = arestfcn(x, armodel_order);

        Pxxar{j} = Pxx_;
        far{j} = f_;
        arcoeff(:,j) = arcoeff_;

        xestar_ = filter(1, arcoeff_, whitenoise);
        xestar_ = xestar_/std(xestar_); % FIXME: forced normalization of AR output. Ideally the noise gain should be known from fundamental relations.
        xestar{j} = xestar_;
        Pxxestar{j} = pwelch(xestar_, window, noverlap, nptswdw, fs, 'onesided');
    end

    % ============
    % Plot results
    % ============
    colors = lines(length(Pxxar));
    
    figure;
    plot(t, x, 'k');
    hold on
    for j=1:length(xestar)
        semilogy(t, xestar{j},'Color', colors(j, :))
    end
    xlabel('Time [s]')
    ylabel(data_unit)
    title(['Sensor data - dataset #' num2str(i)])
    legend([{'Signal (original)'}, arlegendset])
    grid on

    
    figure;
    ax(1) = subplot(211);

    semilogy(f, Pxx, 'k')
    hold on
    for j=1:length(Pxxar)
        semilogy(far{j}, Pxxar{j},'Color', colors(j, :), 'LineWidth', 2)
    end
    xlabel('Frequency [Hz]')
    ylabel(['PSD [(' data_unit ')^2/Hz]'])
    title(['Power spectral density (AR model estimation method) - dataset #' num2str(i)])
    legend([{'Signal (original)'}, arlegendset])
    grid on

    ax(2) = subplot(212);
    semilogy(f, Pxx, 'k')
    hold on
    for j=1:length(Pxxestar)
        semilogy(f, Pxxestar{j},'Color', colors(j, :))
    end
    xlabel('Frequency [Hz]')
    ylabel(['PSD [(' data_unit ')^2/Hz]'])
    title('Power spectral density obtained from filtering white noise with estimated AR models')
    grid on

    linkaxes(ax)

%     figure;
%     plot(fsel, [intrms reveresedintrms], 'LineWidth', 2)
%     grid on
%     xlabel('Frequency [Hz]')
%     ylabel(['Integrated RMS [' data_unit ']'])
%     legend('Direct integration', 'Reverse integration')
    
end