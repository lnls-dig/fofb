function gen_high_level_params_mat(A, ps_pi_fpga_gains_fpath, params_out_fn)
% GEN_HIGH_LEVEL_PARAMS_MAT Generates the high-level parameters matrix file.
%
% A function for generating the high-level parameters matrix file following the
% High-Level FOFB IOC conventions. This file can be directly inputed into the
% High-Level FOFB GUI.
%
% INPUTS:
%   A:                      Open loop transfer functions corrector by corrector
%                           (cell array of dynamical system objects).
%   ps_pi_fpga_gains_fpath: filepath to the power supplies' PI gains obtained by
%                           sysid/ps_pi_tune
%   params_out_fn:          output parameters filename

  % Constants
  Ts = A{1}.Ts;
  ncorr = length(A);
  nbiquads = 4;
  excluded_corr = [1, 80, 81, 160];
  actuator_bw = 10000;

  max_fpga_coeff = 2;

  addpath('sysid');

  ps_pi_fpga_gains = readmatrix(ps_pi_fpga_gains_fpath);

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
  eq_zpk = eqfilt(A,-1,1,-20,1e4,inf,5);
  for i = 1:ncorr
    filters{i} = filter;

    % We don't have fitted models for the excluded correctors
    if ~ismember(i, excluded_corr)
      % Biquad 2-4: Notch filter @ FOFB/2 + Equalization filter

      [eq_i_sos.sos,eq_i_sos.g] = zp2sos(eq_zpk{i}.z{1},eq_zpk{i}.p{1},...
                                         eq_zpk{i}.k);

      % Combine the notch filter at FOFB/2 (1st-order) and the equalization
      % filter (max 5th-order) into the three remaining biquads.
      [b,a] = sos2tf([notch_FOFB_2.sos; eq_i_sos.sos],...
                     notch_FOFB_2.g*eq_i_sos.g);
      sys = tf(b,a,Ts);
      sysr = minreal(sys);
      assert(order(sysr) <= 6);
      assert(isstable(sysr));

      % Convert back to SOS
      [sos,g] = tf2sos(sysr.Numerator{1},sysr.Denominator{1});
      [sos,g] = scale_coeffs_to_fpga(sos,g,max_fpga_coeff);

      filters{i}.sos(2:4, :) = sos;
      filters{i}.g = filters{i}.g*g;
    else
      % Biquad 2: Notch filter @ FOFB/2
      filters{i}.sos(2, :) = notch_FOFB_2.sos;
      filters{i}.g = filters{i}.g*notch_FOFB_2.g;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Plot composed filters' bode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %figure();
  %opts = bodeoptions;
  %opts.FreqUnits = 'Hz';
  %opts.PhaseWrapping = 'on';

  %for i = 1:ncorr
  %  if ~ismember(i, excluded_corr)
  %    [b, a] = sos2tf(filters{i}.sos, filters{i}.g);
  %    sys = tf(b, a, Ts);

  %    % DC gain should be 0 dB
  %    assert((dcgain(sys) - 1.0) < 0.001);

  %    bode(sys, opts);
  %    grid on;
  %    hold on;
  %  end
  %end
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
