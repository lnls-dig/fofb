% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% A function for tuning the power supplies' PI controllers to achieve a
% desired closed-loop bandwidth. This is based on [1].
%
% Credits to Gabriel Ramirez <gabriel.ramirez@lnls.br> for first
% implementing this on 'fofb_ps_sysid.m'.
%
% [1] https://www.ti.com/lit/ta/ssztcg4/ssztcg4.pdf

function ps_pi_fpga_gains = ps_pi_tune(ol_ps_idtf_fpath, Ts, cl_bw)
    ol_ps_idtf = load(ol_ps_idtf_fpath).ol_ps_idtf;
    ncorr = size(ol_ps_idtf, 1);

    figure();

    ps_pi_fpga_gains = zeros(ncorr, 2);
    parfor i = 1:ncorr
        if ~ismember(i, [1, 80, 81, 160])
            fprintf('Corrector #%d... ', i);

            R =  1/dcgain(ol_ps_idtf{i});
            L = -R/pole(ol_ps_idtf{i});
            kp = 2*pi*cl_bw*L;
            ki = 2*pi*cl_bw*R;

            C = tf([kp ki], [1 0]);
            opts = bodeoptions;
            opts.FreqUnits = 'Hz';
            bode(feedback(C*ol_ps_idtf{i}, 1), opts);
            hold on;

            % Hardcoded for now
            quant_gain_v = 2*3.7/(2^16);
            quant_gain_i = 6.25e-5;
            quant_gain = quant_gain_i/quant_gain_v;
            fp_pos = 16;

            ps_pi_fpga_gains(i, :) = [int32(kp*quant_gain*(2^fp_pos)), ...
                                      int32(ki*Ts*quant_gain*(2^fp_pos))];

            fprintf(' done!\n', i);
        end
    end

    ps_pi_fpga_gains_fname = split(ol_ps_idtf_fpath, '-ol_ps_idft.mat');
    ps_pi_fpga_gains_fname = join(ps_pi_fpga_gains_fname, ...
                                  '-ps_pi_fpga_gains.txt');
    writematrix(ps_pi_fpga_gains, ps_pi_fpga_gains_fname{1});
end
