% GEN_HIGH_LEVEL_PARAMS_MAT Generates the high-level parameters matrix file.
%
%   A function for generating the high-level parameters matrix file following
%   the High-Level FOFB IOC conventions. This file can be directly inputed into
%   the High-Level FOFB GUI.
%
%   The arguments are:
%   cl_ps_idtf_fpath: filepath to the power supplies' closed-loop fitted models
%                     obtained by sysid/ps_cl_tfest
%   ps_pi_fpga_gains_fpath: filepath to the power supplies' PI gains file
%                           obtained by sysid/ps_pi_tune
%   params_out_fn: output parameters filename

% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>

function gen_high_level_params_mat(cl_ps_idtf_fpath, ps_pi_fpga_gains_fpath, ...
                                   params_out_fn)
  % Constants
  fs = 48193;
  ncorr = 160;
  nbiquads = 4;
  excluded_corr = [1, 80, 81, 160];
  actuator_bw = 10000;

  addpath('sysid');

  cl_ps_idtf = load(cl_ps_idtf_fpath).cl_ps_idtf;
  ps_pi_fpga_gains = readmatrix(ps_pi_fpga_gains_fpath);

  pass_through_biquad.sos = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
  pass_through_biquad.g = 1.0;

  filter.sos = repmat(pass_through_biquad.sos, nbiquads, 1);
  filter.g = pass_through_biquad.g;

  % Common filters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Biquad 1: Notch filter @ FOFB/4
  notch_FOFB_2 = calc_notch_biquad(0.5, 8);
  filter.sos(1, :) = notch_FOFB_2.sos;
  filter.g = filter.g*notch_FOFB_2.g;

  % Biquad 2: Notch filter @ FOFB/2
  notch_FOFB_4 = calc_notch_biquad(0.9999999999, 8);
  filter.sos(2, :) = notch_FOFB_4.sos;
  filter.g = filter.g*notch_FOFB_4.g;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Specific filters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filters = cell(ncorr, 1);
  for i = 1:ncorr
    filters{i} = filter;

    % We don't have fitted models for the excluded correctors
    if ~ismember(i, excluded_corr)
      cl_ps_bw = bandwidth(cl_ps_idtf{i})/(2*pi);
      % Biquad 3: Pre-emphasis shelf filter
      pre_emph_shelf = calc_shelf_biquad(-cl_ps_bw, -actuator_bw, 1/fs);
      filters{i}.sos(3, :) = pre_emph_shelf.sos;
      filters{i}.g = filters{i}.g*pre_emph_shelf.g;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Plot composed filters' bode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure();
  opts = bodeoptions;
  opts.FreqUnits = 'Hz';
  opts.PhaseWrapping = 'on';

  for i = 1:ncorr
    if ~ismember(i, excluded_corr)
      [b, a] = sos2tf(filters{i}.sos, filters{i}.g);
      sys = tf(b, a, 1/fs);

      % DC gain should be 0 dB
      assert((dcgain(sys) - 1.0) < 0.001);

      bode(sys, opts);
      hold on;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Build high-level parameters matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nparams = 2 + 1 + 5*nbiquads; % 2 for kp and ki,
                                % 1 for filters' accumulated gain and
                                % 5*nbiquads for filters' coefficients.
  params = zeros(ncorr, nparams);

  % Columns 1, 2: kp, ki
  params(:, 1:2) = ps_pi_fpga_gains;
  for i = 1:ncorr
    filter = filters{i};
    filter.sos(:, 4) = [];  % a0 is assumed to be 1
    % Column 3: filters' accumulated gain
    params(i, 3) = filter.g;
    % Columns 4 to (4 + 5*nbiquads) - 1: filters' coefficients
    params(i, 4:end) = reshape(filter.sos', 1, 5*nbiquads);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  writematrix(params, params_out_fn);
end
