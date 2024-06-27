% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% A function for adjusting the power supplies' PI controllers by multiplying its
% gains by the ratio 'desired/measured' bandwidths. This is used as a refinement
% of ps_pi_tune tuning.

function ps_pi_fpga_gains = ps_pi_adjust(cl_ps_idtf_fpath, ...
                                         ps_pi_fpga_gains_fpath, cl_bw)
    cl_ps_idtf = load(cl_ps_idtf_fpath);
    ps_pi_fpga_gains = readmatrix(ps_pi_fpga_gains_fpath);
    ncorr = size(cl_ps_idtf, 1);

    parfor i = 1:ncorr
        if ~ismember(i, [1, 80, 81, 160])
            fprintf('Corrector #%d... ', i);

            factor = (bandwidth(cl_ps_idtf{i})/(2*pi))/cl_bw;
            ps_pi_fpga_gains(i, :) = int32(ps_pi_fpga_gains(i, :)/factor);

            fprintf(' done!\n', i);
        end
    end

    ps_pi_fpga_gains_fname = split(ps_pi_fpga_gains_fpath, ...
                                   '-ps_pi_fpga_gains.txt');
    ps_pi_fpga_gains_fname = join(ps_pi_fpga_gains_fname, ...
                                  '-adj-ps_pi_fpga_gains.txt');
    writematrix(ps_pi_fpga_gains, ps_pi_fpga_gains_fname{1});
end
