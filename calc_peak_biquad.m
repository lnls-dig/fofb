% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% Function for calculating a biquad peak filter given the peak frequency (n0) in
% units of 'fs/2' (range (0, 1)) and the quality factor (Q).

function biquad = calc_peak_biquad(n0, Q)
    [b, a] = iirpeak(n0, n0/Q);
    sysr = minreal(tf(b, a, 1));
    figure();
    opts = bodeoptions;
    opts.PhaseWrapping = 'on';
    bode(sysr, opts);
    figure();
    zplane(sysr.Numerator{1}, sysr.Denominator{1});

    % Convert to SOS
    [biquad.sos, biquad.g] = tf2sos(sysr.Numerator{1},sysr.Denominator{1});

    % Scale the numerator coefficients to better use FPGA's resolution
    % NOTE: The denominator can't be changed because gateware assumes a0 = 1.
    max_fpga_coeff = 2;
    factor = max_fpga_coeff/max(abs(biquad.sos(1:3)));
    biquad.sos(1:3) = factor*biquad.sos(1:3);
    biquad.g = biquad.g/factor;
end
