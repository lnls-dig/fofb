% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>
% Modified by: Guilherme Ricioli  <guilherme.ricioli@lnls.br>

clc; clear;
addpath 'machine/SIRIUS/'

respmat_fpath = 'respmat.mat';
sysid_res_fpath = 'sysid_res.mat';

%% Simulation parameters

% PRBS
prbs_amplitude = 4000;
prbs_lfsr_len = 9;
prbs_step_duration = 4;
prbs_period = (2^prbs_lfsr_len - 1)*prbs_step_duration;
n_prbs_periods = 2;
assert(n_prbs_periods > 1); % 1st period is considered transient

% Orbit
orbit_rms_noise = 300; % in nm

% System identification type ('Sensitivity' or 'Open Loop')
sysid_type = 'Sensitivity';

%% Building FOFB model
fprintf('Building FOFB model...\n');

M = load(respmat_fpath).mat_d';
A = load(sysid_res_fpath).sys;

excluded_corr = [1 80 81 160];

Ts = A{2}.Ts;

tic
for i=1:size(A,1)
    if ismember(i,excluded_corr)
        A{i} = tf(0,1,'Ts',Ts);
    else
        A{i} = idtf(A{i}/dcgain(idtf(A{i})),'Ts',Ts);
    end
end

fofb_type.bpm_sector = 'M1M2C2C3';
fofb_type.corr_sector = 'M1M2C2C3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = excluded_corr;
[bpm_idx,corr_idx] = fofb_idx(fofb_type);

M = M(bpm_idx,corr_idx);
Mc = pinv(M);
K = (1/Ts)*[0.12*ones(size(Mc,1)/2,1); 0.166*ones(size(Mc,1)/2,1)];
A = A(corr_idx);
[P,G,C] = ofbmdl(M,Mc,K,A);

fprintf('Elapsed time: %.2f s\n',toc);

%% Building PRBS input cube
fprintf('Building PRBS input cube...\n');

prbs_frest = frest.PRBS('Amplitude',prbs_amplitude,...
                   'Ts',Ts,...
                   'Order',prbs_lfsr_len,...
                   'NumPeriods',n_prbs_periods,...
                   'UseWindow','off');
prbs_ts = generateTimeseries(prbs_frest);
prbs = repelem(prbs_ts.Data,prbs_step_duration);
prbs = prbs - mean(prbs);

T=0:Ts:(length(prbs) - 1)*Ts;

% Defining U, V, u and y
[U,S,V] = svd(M);

n_modes = min(size(U,1),size(V,1));
n = size(U,1);
y = zeros(length(T),n,n_modes);
if strcmp(sysid_type,'Sensitivity')
    sys = P('yd',{'d','n'});
    m = n;
    u = zeros(length(T),m,n_modes);
    for i=1:n_modes
        u(:,:,i) = U(:,i)'.*prbs;
    end
elseif strcmp(sysid_type,'Open Loop')
    sys = G;
    m = size(V,1);
    u = zeros(length(T),m,n_modes);
    for i=1:n_modes
        u(:,:,i) = V(:,i)'.*prbs;
    end
end

fprintf('Elapsed time: %.2f s\n',toc);

%% Simulating system
fprintf('Simulating system...\n');

freqs = rfftfreq(prbs_period, Ts);

Yu = ones(length(freqs),m,n_modes);
Yy = ones(length(freqs),n,n_modes);

snr = ones(m,n_modes);

% Iterates for each excited mode
for k=1:n_modes
    fprintf('Mode: %d\n',k);

    if strcmp(sysid_type,'Sensitivity')
        noise = orbit_rms_noise.*randn(size(y,1),size(y,2));
        y(:,:,k) = lsim(sys,[u(:,:,k) noise],T);
    elseif strcmp(sysid_type,'Open Loop')
        y(:,:,k) = lsim(sys,u(:,:,k),T);
        noise = orbit_rms_noise.*randn(size(y,1),1);
        y(:,i,k) = y(:,i,k) + noise;
    end

    % Iterates in each input/output
    for i=1:n
        % Building Yu arrays
        if(i<=m)
            u_avg = reshape(u(:,i,k),prbs_period,n_prbs_periods);
            % 1st PRBS period is considered transient
            u_avg = mean(u_avg(:,2:end),2);
            Yu(:,i,k) = rfft(u_avg);
        end

        % Building Yy arrays
        y_avg = reshape(y(:,i,k),prbs_period,n_prbs_periods);
        % 1st PRBS period is considered transient
        y_avg = mean(y_avg(:,2:end),2);
        Yy(:,i,k) = rfft(y_avg);

        % Calculating SNR
        % if strcmp(sysid_type,'Open Loop')
        %    snr(i,k) = 20*log10(rms(y(:,i,k))/rms(noise));
        % elseif strcmp(sysid_type,'Sensitivity')
        %    snr(i,k) = 20*log10(rms(y(:,i,k))/rms(noise(:,i)));
        % end
    end
end

fprintf('Elapsed time: %.2f s\n',toc);

% Clear unused variables to free up space
clear U V u y prbs T sys;

%% Calculating Response Matrices and their SVDs
fprintf('Calculating Response Matrices and their SVDs...\n');

respmat = zeros(n,m,length(freqs));
exp_sigma = zeros(n,m,length(freqs));

for f=1:length(freqs)
    respmat(:,:,f) = squeeze(Yy(f,:,:))*pinv(squeeze(Yu(f,:,:)));
    [~,exp_sigma(:,:,f),~] = svd(squeeze(respmat(:,:,f)));
end

fprintf('Elapsed time: %.2f s\n',toc);

%% Plotting results
fprintf('Plotting results...\n');

if strcmp(sysid_type,'Sensitivity')
    sys = P('yd','d');
elseif strcmp(sysid_type,'Open Loop')
    sys = G;
end

figure;
for k=1:n_modes
    % We removed DC from PRBS and orbit signals, so start from frequency
    % index 2
    semilogx(freqs(2:end), ...
             20*log10(squeeze(exp_sigma(k,k,2:end))),'Color','#D95319');
    hold on;
end
legend({'Simulated'})
h = sigmaplot(sys);
plotoptions = sigmaoptions;
plotoptions.FreqUnits = 'Hz';
plotoptions.Grid = 'on';
setoptions(h,plotoptions);

fprintf('Elapsed time: %.2f s\n',toc);

%% Savings results
fprintf('Savings results...\n');

save(['exp_sigma_sim-', sysid_type], 'freqs', 'respmat', 'exp_sigma');

fprintf('Elapsed time: %.2f s\n',toc);
fprintf('done!\n');

%% Function definitions
function f = rfftfreq(L, Ts)
    Fs = 1/Ts;
    n = idivide(int32(L),int32(2),'ceil');
    f = (Fs/L)*(0:double(n) - 1);
end

function Y = rfft(X)
    L = length(X);
    n = idivide(int32(L),int32(2),'ceil');
    Y = fft(X)/L;
    Y = Y(1:n);
end
