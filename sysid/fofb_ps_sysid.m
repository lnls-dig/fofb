% This script was created to analyse the data acquired from the power supplies 
% in the machine shifts conduced throughout may 2023 as part of the system
% identification process for the FOFB system. 

% It aims at fullfilling the following goals:
%
% 1) Identify the open loop response of each one of the power supplies that
%    are part of the FOFB system and estimate the inductances and
%    resistances seem by each one of them;
%
% 2) Based on the open loop identification, calculate the PI gains that
%    should equalize the responses of each of the power supplies while
%    providing the desired closed loop bandwidth for them;
%
% 3) Once the calculated/equalizing PI gains are applied to the power
%    supplies' controllers, identify the closed loop responses and verify if
%    the estimated models from the closed loop data agree with the expected
%    responses.
%
% For the purposes of identification and equalization of the power
% supplies, step signals were used as probing signals


%% Defines file locations and names for the open loop acquisitions

open_loop_base_path = '~/Machine_Studies/2023-05-29-SI_FOFB_SYSID/shared-mount/mat/voltage-step-300mv/';
open_loop_acqs_fname = 'lamp-300mv-01.mat';
open_loop_acqs = load([open_loop_base_path open_loop_acqs_fname]);

%% Extract input and output data for open loop estimation

open_loop_out = open_loop_acqs.current;
open_loop_in = open_loop_acqs.voltage;

%% Definitions of some general parameters for the script

ncorr = size(open_loop_out, 2);
Ts = 1.1556e-6;
num_histogram_bins = 20;

% Used mainly in the plotting sections of the script
hor_corr_idxs = 1:80;
vert_corr_idxs = 81:160;
excluded_corrs = [1, 80, 81, 160];

%% Model parameters definitions for the open loop estimation

np_ol = 1;
nz_ol = 0;
dly_ol = 2*Ts;
idx_crop = 1:70e3;

%% Open loop system identification of the power supplies

ol_estimated_power_supplies = cell(ncorr, 1);
open_loop_bws = zeros(ncorr,1);

for i=1:ncorr
    tic;
    ol_estimated_power_supplies{i} = tfest(iddata(open_loop_out(idx_crop,i),...
                                  open_loop_in(idx_crop,i), Ts), np_ol, nz_ol, dly_ol);
    fprintf('Corrector #%d...', i);
    elapsed_time = toc;
    fprintf(' %0.3f s.\n', elapsed_time);
    
    open_loop_bws(i) = bandwidth(ol_estimated_power_supplies{i})/(2*pi); % [Hz]
end

%% Plots for the fit percentages of the estimated models
ol_fit_percentage = zeros(ncorr, 1);

for i = 1:ncorr
   ol_fit_percentage(i) = ol_estimated_power_supplies{i}.Report.Fit.FitPercent; 
end 

ol_fit_percentage(excluded_corrs) = NaN;

figure;
subplot(211);
plot(ol_fit_percentage(hor_corr_idxs));
ylabel('Fit percentage (%)')
hold on;
title('Horizontal')

subplot(212);
plot(ol_fit_percentage(vert_corr_idxs));
hold on;
ylabel('Fit percentage (%)')
xlabel('Corrector Index')
title('Vertical')

sgtitle('Open loop estimation fit percentage')

%% Plots of the estimated bandwidths of the open loop power supplies

open_loop_bws(excluded_corrs) = NaN;

% ------- Values -------
figure;
subplot(211);
plot(open_loop_bws(hor_corr_idxs));
hold on;
title('open loop FCHs Bandwidths [Hz]')

subplot(212);
plot(open_loop_bws(vert_corr_idxs));
hold on;
title('open loop FCVs Bandwidths [Hz]')

% ------- Histograms -------
figure;
subplot(211);
histogram(open_loop_bws(hor_corr_idxs));
hold on;
title('open loop FCHs Bandwidths [Hz]')

subplot(212);
histogram(open_loop_bws(vert_corr_idxs));
hold on;
title('open loop FCVs Bandwidths [Hz]')

%% Estimation of the inductances and resistances of the correctors

R = zeros(ncorr,1);
L = zeros(ncorr,1);

for i=1:ncorr
    L(i) = 1/ol_estimated_power_supplies{i}.Numerator;
    R(i) = -L(i)*pole(ol_estimated_power_supplies{i});
end

%% Plots of the estimated resistances and inductances

L_aux = L;
R_aux = R;

L_aux(excluded_corrs) = NaN;
R_aux(excluded_corrs) = NaN;

% ------------- Horizontal ---------------
figure;
subplot(211);
plot(L_aux(hor_corr_idxs)/1e-3);
hold on;
ylabel('Inductance [mH]')
title('Inductances')

subplot(212);
plot(R_aux(hor_corr_idxs));
hold on;
ylabel('Resistance [Ohm]')
xlabel('Corrector Index')
title('Resistances')

sgtitle('Horizontal')

% ------------- Vertical ---------------
figure;
subplot(211);
plot(L_aux(vert_corr_idxs)/1e-3);
hold on;
ylabel('Inductance [mH]')
title('Inductances')

subplot(212);
plot(R_aux(vert_corr_idxs));
ylabel('Resistance [Ohm]')
xlabel('Corrector Index')
title('Resistances')
hold on;
title('Resistances')

sgtitle('Vertical')
%% Histograms of the estimaed resistances and inductances

L_aux = L;
R_aux = R;

L_aux(excluded_corrs) = NaN;
R_aux(excluded_corrs) = NaN;

% ------------- Horizontal ---------------
figure; 
subplot(211); 
histogram(L_aux(hor_corr_idxs)/1e-3, num_histogram_bins);
title('Inductances [mH]')

subplot(212);
histogram(R_aux(hor_corr_idxs), num_histogram_bins);
title('Resistences [Ohm]')

sgtitle('Horizontal')

% ------------- Vertical ---------------
figure; 
subplot(211); 
histogram(L_aux(vert_corr_idxs)/1e-3, num_histogram_bins);
title('Inductances [mH]')

subplot(212);
histogram(R_aux(vert_corr_idxs), num_histogram_bins);
title('Resistences [Ohm]')

sgtitle('Vertical')

%% Design of the closed loop PI gains

clbw = 4e+3; %(Hz)

Kp = zeros(length(L), 1);
Ki = zeros(length(R), 1);

for i=1:length(L)
    Kp(i) = 2*pi*clbw*L(i);
    Ki(i) = 2*pi*clbw*R(i);
end

%% Calculates the theoretical/expected closed loop responses and their bandwidths

sim_current_closed_loop = {};
sim_voltage_actuation = {};
sim_closed_loop_bws = zeros(ncorr, 1);

for i=1:length(Kp)
    pi_controller = pid(Kp(i), Ki(i));
    sim_current_closed_loop{i} = feedback(pi_controller*ol_estimated_power_supplies{i}, 1);
    sim_voltage_actuation{i} = feedback(pi_controller, ol_estimated_power_supplies{i});
    sim_closed_loop_bws(i) = bandwidth(sim_current_closed_loop{i})/(2*pi);
end

%% Plots for the expected closed loop bandwidths

sim_closed_loop_bws(excluded_corrs) = NaN;

figure
subplot(211);
plot(sim_closed_loop_bws(hor_corr_idxs)/1e3);
hold on;
ylabel('Bandwidth [KHz]')
title('Expected FCHs Bandwidths')

subplot(212);
plot(sim_closed_loop_bws(vert_corr_idxs)/1e3);
hold on;
ylabel('Bandwidth [KHz]')
xlabel('Corrector Index');
title('Expected FCVs Bandwidths')

%% Defines closed loop acquisitions file locations and names

closed_loop_base_path = '~/Machine_Studies/2023-05-30-SI_FOFB_SYSID/shared-mount/mat/current-step-30ma-5khz/';

%% Extract input and output data for closed loop estimation
num_cl_acqs = 200;
num_pts = 3000;
closed_loop_out_avg = zeros(num_pts, ncorr);

for i=1:num_cl_acqs    
    
    closed_loop_acqs_fname = sprintf('lamp-%03d.mat', i);
    closed_loop_acqs = load([closed_loop_base_path closed_loop_acqs_fname]);
    closed_loop_out_avg =  closed_loop_out_avg + closed_loop_acqs.current/num_cl_acqs;
end

%% Closed loop input definition

% Since we are not able to acquire the probing input signal at the
% present moment when we are exciting the closed loop with a current
% signal. We generate it artificially based on the supposed current step
% that should have been aplied

step_start = 101;
step_amplitude = 30e-3;
signal_duration = 3000;

current_step_input = zeros(signal_duration, 1);
current_step_input(step_start: end) = step_amplitude;

%% Model parameters definitions for the open loop estimation

np_cl = 2;
nz_cl = 1;
dly_cl = 2*Ts;

%% Closed loop system identification of the power supplies

meas_closed_loop_bws = zeros(ncorr,1);
cl_estimated_power_supplies = cell(ncorr, 1);

for i=1:ncorr
    tic;
    cl_estimated_power_supplies{i} = tfest(iddata(closed_loop_out_avg(:,i), ... 
                                     current_step_input, Ts), np_cl, nz_cl, dly_cl);
    fprintf('Corrector #%d...', i);
    elapsed_time = toc;
    fprintf(' %0.3f s.\n', elapsed_time);
    
    meas_closed_loop_bws(i) = bandwidth(cl_estimated_power_supplies{i})/(2*pi); % [Hz]
end

%% Plots for the fit percentages of the closed loop estimated models
cl_fit_percentage = zeros(ncorr, 1);

for i = 1:ncorr
   cl_fit_percentage(i) = cl_estimated_power_supplies{i}.Report.Fit.FitPercent; 
end 

cl_fit_percentage(excluded_corrs) = NaN;

figure;
subplot(211);
plot(cl_fit_percentage(hor_corr_idxs));
ylabel('Fit percentage (%)')
hold on;
title('Horizontal')

subplot(212);
plot(cl_fit_percentage(vert_corr_idxs));
hold on;
ylabel('Fit percentage (%)')
xlabel('Corrector Index')
title('Vertical')

sgtitle('Closed loop estimation fit percentage')

%% Plots of the estimated bandwidths of the closed loop power supplies

meas_closed_loop_bws(excluded_corrs) = NaN;

% ------ Values ------
figure;
subplot(211);
plot(meas_closed_loop_bws(hor_corr_idxs)/1e3);
hold on;
ylabel('Bandwidth [KHz]')
title('closed loop FCHs Bandwidths')

subplot(212);
plot(meas_closed_loop_bws(vert_corr_idxs)/1e3);
hold on;
ylabel('Bandwidth [KHz]')
xlabel('Corrector Index')
title('closed loop FCVs Bandwidths [kHz]')

% ------ Histograms ------
figure;
subplot(211);
histogram(meas_closed_loop_bws(hor_corr_idxs)/1e3, num_histogram_bins);
hold on;
title('closed loop FCHs Bandwidths [kHz]')

subplot(212);
histogram(meas_closed_loop_bws(vert_corr_idxs)/1e3, num_histogram_bins);
hold on;
title('closed loop FCVs Bandwidths [kHz]')

%% Convert pi gains to fpga gains

voltage_count = 2*3.7/(2^16);
current_count = 6.25e-5;
conversion_factor = voltage_count/current_count;

Kp_fpga = zeros(length(Kp), 1);
Ki_fpga = zeros(length(Ki), 1);

for i=1:length(Kp)
    Kp_fpga(i) = (Kp(i)*(2^16))/conversion_factor;
    Ki_fpga(i) = (Ki(i)*Ts*(2^16))/conversion_factor;
end

%% 
kp = 5000000;
ki = 2000;
Vref = 3.7;
curr_Ts = 1.15e-6;
% curr_Ts = Ts;

conversion_factor = (2*Vref/2^16)/6.25e-5;

kp_physical = (kp*conversion_factor)/(2^16);
ki_physical = (ki*conversion_factor)/(curr_Ts*(2^16));