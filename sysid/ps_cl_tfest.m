% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>
%
% A function for estimating closed-loop models FOFB power supplies.
%
% Credits to Gabriel Ramirez <gabriel.ramirez@lnls.br> for first
% implementing this on 'fofb_ps_sysid.m'.

function cl_ps_idtf = ps_cl_tfest(cl_acq_fpath, slice_to_fit, Ts, np, nz, ...
                                  iodelay)
    cl_acq = load(cl_acq_fpath);
    ncorr = size(cl_acq.current, 2);

    figure();

    % Hardcoded for now
    r = zeros(length(slice_to_fit), 1);
    r(101:end) = 0.02;

    cl_ps_iddata = cell(ncorr, 1);
    cl_ps_idtf = cell(ncorr, 1);
    parfor i = 1:ncorr
        if ~ismember(i, [1, 80, 81, 160])
            fprintf('Corrector #%d... ', i);

            current = smoothdata(cl_acq.current(slice_to_fit, i), 'movmean', ...
                                 200);
            cl_ps_iddata{i} = iddata(current, r, Ts);
            tic;
            cl_ps_idtf{i} = tfest(cl_ps_iddata{i}, np, nz, iodelay);
            fprintf('took %0.3fs (fit: %0.2f).\n', toc, ...
                    cl_ps_idtf{i}.Report.Fit.FitPercent);

            % figure();
            % t = (slice_to_fit - 1)*Ts;
            % plot(t, cl_acq.current(slice_to_fit, i), ...
            %      t, lsim(cl_ps_idtf{i}, r, t));
            % legend({'real', 'model'})
            %
            % figure();
            % opts = bodeoptions;
            % opts.FreqUnits = 'Hz';
            % bode(cl_ps_idtf{i}, opts);

            opts = bodeoptions;
            opts.FreqUnits = 'Hz';
            bode(cl_ps_idtf{i}, opts);
            hold on;
        end
    end

    cl_ps_idft_fname = split(cl_acq_fpath, '.mat');
    cl_ps_idft_fname = join(cl_ps_idft_fname, '-cl_ps_idft.mat');
    save(cl_ps_idft_fname{1}, 'cl_ps_idtf');
end
