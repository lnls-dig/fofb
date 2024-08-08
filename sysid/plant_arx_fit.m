% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

function [sys,fit] = plant_arx_fit(file, arx_parameters, orbit_plane, acquisition_period)
% plant_arx_fit
%
% Provides an ARX fit to the studied system with PRBS excitation
%
% [sys,fit] = plant_arx_fit(file, arx_parameters, orbit_plane, acquisition_period)
%
% INPUTS:
%   file:               File in format .mat with input and output data from
%                       sysid experiment after processing it with h5tomat.
%   arx_parameters:     Parameters for the order of the ARX estimation.
%                       Expected in the format [Na Nb Ts], where Na and Nb
%                       are the order of polynomials A and B respectively
%                       and Ts is the number of delay samples considered.
%                       More details may be found in the documentation for
%                       the arx function.
%   orbit_plane:        Orbit plane to fit the system. Available options:
%                       'X' or 'Y'.
%   acquisition_period: (Optional input) If the studied aquisition time is
%                       intended to be less than what the arrays contain,
%                       this variable may be set to limit it.
%
% OUTPUTS:
%   sys:                Resulting system in idpoly format.
%   fit:                Fit percent of obtained system.

%Loads .mat files
A = load(file,'data');
input_signal = A.data.input_signal;
if orbit_plane=='X'
    orbit_data = double(A.data.orbx);
elseif orbit_plane=='Y'
    orbit_data = double(A.data.orby);
else
    print('Error: orbit_plane value not allowed')
end
sampling_frequency = A.data.sampling_frequency;
prbs_lfsr_len = A.data.prbs_lfsr_len(1);
prbs_step_duration = A.data.prbs_step_duration(1);
prbs_mov_avg_taps = double(A.data.prbs_mov_avg_taps);
%Finds period length for PRBS excitation
n=(2^(prbs_lfsr_len)-1)*prbs_step_duration;

%Removes transient (removes 4 periods)
transient_samples = 4*n;
assert(length(input_signal)>transient_samples,"Could not remove transient: input_signal array is too small");
input_signal = input_signal(transient_samples:end);
assert(length(orbit_data)>transient_samples,"Could not remove transient: orbit_data array is too small");
orbit_data = orbit_data(transient_samples:end);

%Creates two disjoint datasets: estimation and validation data
prbs_periods = length(orbit_data)/n;
h=int64(prbs_periods/2);
input_signal = input_signal(1:h*n);
orbit_data = orbit_data(1:h*n);

%Reduction in aquisition period (optional parameter)
T_aq=single(h)/sampling_frequency;
if exist('aquisition_period', 'var')
    k=floor(acquisition_period*h);
    fprintf('Total time considered: %f s\n', T_aq*single(k)/single(h))
    input_signal = input_signal(1:k);
    orbit_data = orbit_data(1:k);
end

%Applied filters
fir_moving_average = dfilt.dffir(ones(1,2^prbs_mov_avg_taps)/2^prbs_mov_avg_taps);
%Low pass filter: f3db = 231.4866 mHz. Obtained from measure(iir_lowpass)
%iir_lowpass = design(fdesign.lowpass(0.23, 0.25, 1, 80, 1), 'butter', 'MatchExactly', 'stopband');
input_signal = filter(fir_moving_average, input_signal);
orbit_data = filter(fir_moving_average, orbit_data);
%currdata = filter(iir_lowpass, currdata);
%orbit_data = filter(iir_lowpass, orbit_data);

%Removes average
input_signal=input_signal-mean(input_signal);
orbit_data=orbit_data-mean(orbit_data);

%Adjusts size
%currdata = currdata(1:floor(size(currdata,1)/n)*n, :);
%orbit_data = orbit_data(1:floor(size(orbit_data,1)/n)*n, :);

%Alternative solution code snippet
%Defines smaller arrays
%sublengths = ones(1, floor(numel(currdata)/n))*n;
%currdata_div = mat2cell(currdata', 1, sublengths);
%aux=zeros(length(sublengths),sublengths(1));
%for i=1:numel(sublengths)
%    aux(i,:)=currdata_div{1,i};
%    currdata_avg=mean(aux,1)';
%end
%sublengths = ones(1, floor(numel(orbit_data)/n))*n;
%orbit_data_div = mat2cell(orbit_data', 1, sublengths);
%for i=1:numel(sublengths)
%    aux(i,:)=orbit_data_div{1,i};
%    orbit_data_avg=mean(aux,1)';
%end

orbit_data_div = reshape(orbit_data,n,[]);
orbit_data_avg = mean(orbit_data_div,2);

input_signal_div = reshape(input_signal,n,[]);
input_signal_avg = mean(input_signal_div,2);

%ARX fit
sys = arx(input_signal_avg, orbit_data_avg, arx_parameters, arxOptions('Focus', 'Simulation', 'EnforceStability', true));
fit=sys.Report.Fit.FitPercent;