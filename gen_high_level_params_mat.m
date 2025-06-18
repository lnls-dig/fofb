% GEN_HIGH_LEVEL_PARAMS_MAT Generates the high-level parameters matrix file.
%
%   A function for generating the high-level parameters matrix file following
%   the High-Level FOFB IOC conventions. This file can be directly inputed into
%   the High-Level FOFB GUI.
%
%   The arguments are:
%   cl_ps_idtf_fpath: filepath to the power supplies' closed-loop fitted models
%                     obtained by sysid/ps_cl_tfest
%   ps_pi_fpga_gains_fpath: filepath to the power supplies' PI gains obtained by
%                           sysid/ps_pi_tune
%   F_fpath: filepath to the cell array of dimensions NCORR x 1 of objects 'tf'
%            representing the filters to be incorporated
%   params_out_fn: output parameters filename

% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Guilherme Ricioli <guilherme.ricioli@lnls.br>

function gen_high_level_params_mat(cl_ps_idtf_fpath, ps_pi_fpga_gains_fpath, ...
                                   F_fpath, params_out_fn)
  % Constants
  fs = 48193;
  ncorr = 160;
  nbiquads = 4;
  excluded_corr = [1, 80, 81, 160];
  actuator_bw = 10000;

  addpath('sysid');

  cl_ps_idtf = load(cl_ps_idtf_fpath).cl_ps_idtf;
  ps_pi_fpga_gains = readmatrix(ps_pi_fpga_gains_fpath);
  F = load(F_fpath).F;

  pass_through_biquad.sos = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
  pass_through_biquad.g = 1.0;

  filter.sos = repmat(pass_through_biquad.sos, nbiquads, 1);
  filter.g = pass_through_biquad.g;

  % Common biquads
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Biquad 1: Notch filter @ FOFB/4
  notch_FOFB_4 = calc_notch_biquad(0.5, 8);
  filter.sos(1, :) = notch_FOFB_4.sos;
  filter.g = filter.g*notch_FOFB_4.g;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Specific biquads
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filters = cell(ncorr, 1);

  notch_FOFB_2 = calc_notch_biquad(0.9999999999, 8);
  for i = 1:ncorr
    filters{i} = filter;

    % We don't have fitted models for the excluded correctors
    if ~ismember(i, excluded_corr)
      cl_ps_bw = bandwidth(cl_ps_idtf{i})/(2*pi);

      % Biquad 2: Notch filter @ FOFB/2 + Pre-emphasis shelf filter
      % Since both the notch filter at FOFB/2 and the pre-emphasis shelf filter
      % are first-order systems, they can be combined into a single biquad.

      pre_emph_shelf = calc_shelf_biquad(-cl_ps_bw, -actuator_bw, 1/fs);
      [b,a] = sos2tf([notch_FOFB_2.sos; pre_emph_shelf.sos], ...
                     notch_FOFB_2.g*pre_emph_shelf.g);
      sys = tf(b,a,1/fs);
      sysr = minreal(sys);
      assert(order(sysr) == 2);
      assert(isstable(sysr));

      % Convert to SOS
      [biquad.sos, biquad.g] = tf2sos(sysr.Numerator{1},sysr.Denominator{1});

      % Scale the numerator coefficients to better use FPGA's resolution
      % NOTE: The denominator can't be changed because gateware assumes a0 = 1.
      max_fpga_coeff = 2;
      factor = max_fpga_coeff/max(abs(biquad.sos(1:3)));
      biquad.sos(1:3) = factor*biquad.sos(1:3);
      biquad.g = biquad.g/factor;

      filters{i}.sos(2, :) = biquad.sos;
      filters{i}.g = filter.g*biquad.g;

      % Biquads 3 and 4: Equalization filter
      assert(order(F{i}) <= 4);
      assert(isstable(F{i}));

      % Convert to SOS
      [biquad.sos, biquad.g] = tf2sos(F{i}.Numerator{1},F{i}.Denominator{1});

      % Scale the numerator coefficients to better use FPGA's resolution
      % NOTE: The denominator can't be changed because gateware assumes a0 = 1.
      max_fpga_coeff = 2;
      for j = 1:size(biquad.sos,1)
        factor = max_fpga_coeff/max(abs(biquad.sos(j,1:3)));
        biquad.sos(j,1:3) = factor*biquad.sos(j,1:3);
        biquad.g = biquad.g/factor;

        filters{i}.sos(2+j,:) = biquad.sos(j,:);
      end

      filters{i}.g = filters{i}.g*biquad.g;
    else
      % Biquad 2: Notch filter @ FOFB/2
      filters{i}.sos(2, :) = notch_FOFB_2.sos;
      filters{i}.g = filters{i}.g*notch_FOFB_2.g;
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
      grid on;
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
