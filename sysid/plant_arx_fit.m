% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>
% Modified by: Guilherme Ricioli <guilherme.ricioli@lnls.br>

function sys = plant_arx_fit(fpath, arx_params, n_prbs_T_to_use)
% plant_arx_fit
%
% Provides an ARX fit to the studied system with PRBS excitation
%
% sys = plant_arx_fit(fpath, arx_params, n_prbs_T_to_use)
%
% INPUTS:
%   fpath:              Filepath to the plant PRBS acquistion MATLAB object
%                       generated by prbs_ol_acq_h5_to_mat.m script
%   arx_params:         ARX estimation parameters in the format [na nb nk] (see
%                       'arx' documentation)
%   n_prbs_T_to_use:    (optional parameter) Number of PRBS periods to be
%                       considered
%
% OUTPUTS:
%   sys:                Resulting system in idpoly format

% Loads .mat fpath
prbs_ol_acq = load(fpath);
Ts = 1/prbs_ol_acq.data.sampling_frequency;
prbs_u = prbs_ol_acq.data.prbs_data;
prbs_lfsr_len = prbs_ol_acq.data.prbs_lfsr_len(1);
prbs_step_duration = prbs_ol_acq.data.prbs_step_duration(1);
prbs_mov_avg_taps = double(prbs_ol_acq.data.prbs_mov_avg_taps);
bpm_y = prbs_ol_acq.data.orb(prbs_ol_acq.data.bpm_idx_max_response, :);

% PRBS period length
prbs_T = (2^(prbs_lfsr_len) - 1)*prbs_step_duration;

% Removes transient
n_prbs_T_transient = 4;
assert(length(prbs_u) > prbs_T*n_prbs_T_transient, ...
       "Could not remove transient: dataset is too small");
prbs_u = prbs_u(prbs_T*n_prbs_T_transient + 1:end);
bpm_y = bpm_y(prbs_T*n_prbs_T_transient + 1:end);

% Guarantees that there's an integer amount of PRBS periods in datasets
rem = mod(length(prbs_u), prbs_T);
prbs_u = prbs_u(1:end - rem);
bpm_y = bpm_y(1:end - rem);

% Reduces dataset (optional parameter)
if exist('n_prbs_T_to_use', 'var')
    n_p = length(prbs_u)/prbs_T;
    assert(n_prbs_T_to_use <= n_p, ...
           "Could not reduce dataset: not enough PRBS periods");
    prbs_u = prbs_u(1:n_prbs_T_to_use*prbs_T);
    bpm_y = bpm_y(1:n_prbs_T_to_use*prbs_T);
end

% Mimics the moving average filter that's applied to the PRBS excitation in
% gateware so we don't end up identifying it
mov_avg_filt = dfilt.dffir(ones(1, 2^prbs_mov_avg_taps)/2^prbs_mov_avg_taps);
prbs_u = filter(mov_avg_filt, prbs_u);

% Averaging
prbs_u_avg = mean(reshape(prbs_u, prbs_T, []), 2);
bpm_y_avg = mean(reshape(bpm_y, prbs_T, []), 2);

% Removes the mean value
prbs_u_avg = prbs_u_avg - mean(prbs_u_avg);
bpm_y_avg = bpm_y_avg - mean(bpm_y_avg);

% ARX fit
plant_iddata = iddata(bpm_y_avg, prbs_u_avg, Ts);
sys = arx(plant_iddata, arx_params, ...
          arxOptions('Focus', 'Simulation', 'EnforceStability', true));
