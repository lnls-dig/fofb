% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% A function for estimating open-loop models for FOFB power supplies.
%
% Credits to Gabriel Ramirez <gabriel.ramirez@lnls.br> for first implementing
% this on 'fofb_ps_sysid.m'.

function ol_ps_idtf = ps_ol_tfest(ol_acq_fpath, slice_to_fit, Ts, np, nz, ...
                                  iodelay)
    ol_acq = load(ol_acq_fpath);
    ncorr = size(ol_acq.voltage, 2);

    figure();

    ol_ps_iddata = cell(ncorr, 1);
    ol_ps_idtf = cell(ncorr, 1);
    parfor i = 1:ncorr
        if ~ismember(i, [1, 80, 81, 160])
            fprintf('Corrector #%d ', i);

            ol_ps_iddata{i} = iddata(ol_acq.current(slice_to_fit, i), ...
                                     ol_acq.voltage(slice_to_fit, i), Ts);
            tic;
            ol_ps_idtf{i} = tfest(ol_ps_iddata{i}, np, nz, iodelay);
            fprintf('took %0.3fs (fit: %0.2f).\n', toc, ...
                    ol_ps_idtf{i}.Report.Fit.FitPercent);

            % figure();
            % t = (slice_to_fit - 1)*Ts;
            % plot(t, ol_acq.current(slice_to_fit, i), ...
            %      t, lsim(ol_ps_idtf{i}, ol_acq.voltage(slice_to_fit, i), t));
            % legend({'real', 'model'})
            % close();
            %
            % figure();
            % opts = bodeoptions;
            % opts.FreqUnits = 'Hz';
            % bode(ol_ps_idtf{i}, opts);

            opts = bodeoptions;
            opts.FreqUnits = 'Hz';
            bode(ol_ps_idtf{i}, opts);
            hold on;
          end
    end

    ol_ps_idft_fname = split(ol_acq_fpath, '.mat');
    ol_ps_idft_fname = join(ol_ps_idft_fname, '-ol_ps_idft.mat');
    save(ol_ps_idft_fname{1}, 'ol_ps_idtf');
end
