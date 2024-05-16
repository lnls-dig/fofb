% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% Function for calculating a biquad shelf filter given the zero (z) and pole (p)
% (both in Hz) of the corresponding analog filter and the sampling time (Ts).
% NOTE: The filter is designed to have DC gain of 0dB.

function biquad = calc_shelf_biquad(z, p, Ts)
    wz = 2*pi*z;
    wp = 2*pi*p;
    k = p/z;
    % Analog filter
    sysc = tf(zpk(wz, wp, k));
    figure();
    opts = bodeoptions;
    opts.FreqUnits = 'Hz';
    opts.PhaseWrapping = 'on';
    bode(sysc, opts);
    hold on;

    % Discretize using Tustin's method
    sysd = c2d(sysc, Ts, 'Tustin');
    bode(sysd, opts);
    legend({'continuous', 'discrete'});
    figure()
    zplane(sysd.Numerator{1}, sysd.Denominator{1});

    % Convert to SOS
    [biquad.sos, biquad.g] = tf2sos(sysd.Numerator{1}, sysd.Denominator{1});

    % Scale the numerator coefficients to better use FPGA's resolution
    % NOTE: The denominator can't be changed because gateware assumes a0 = 1.
    max_fpga_coeff = 2;
    factor = max_fpga_coeff/max(abs(biquad.sos(1:3)));
    biquad.sos(1:3) = factor*biquad.sos(1:3);
    biquad.g = biquad.g/factor;
end
