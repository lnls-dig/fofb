function [orb, orb_envelope, f] = fofb_orb_by_freq(fa_data, f)
% [orb, orb_envelope, f] = fofb_orb_by_freq(fa_data, f)

fft_result = fft(fa_data.bpm_readings);
[npts, n_bpms] = size(fa_data.bpm_readings);

bpm_position_fft_phase = angle(fft_result);
bpm_position_fft_mag = 2/npts*abs(fft_result);

% Calculate sampling frequency and convert to Hz
Fs = 1e9/double(mean(diff(fa_data.time)));
freq = linspace(0, Fs, npts+1);
freq = freq(1:end-1);

orb = zeros(n_bpms, length(f));
orb_envelope = zeros(n_bpms, length(f));
for i=1:length(f)
    mag_by_freq = interp1(freq, bpm_position_fft_mag, f(i),'nearest');

    phase_by_freq = interp1(freq, bpm_position_fft_phase, f(i),'nearest');
    phase_by_freq = phase_by_freq - repmat(phase_by_freq(:,1), 1, size(phase_by_freq,2));

    orb(:,i) = mag_by_freq.*cos(phase_by_freq);
    orb_envelope(:,i) = mag_by_freq;
end