%FOFB model with RF transfer function (rf phase noise -> BPM) implemented
%within ofbmdl_rf.m function

clc;clear;
%FOFB parameters
load respmat.mat
load sysid_res.mat
M = mat_d';
Ts=1/48193;
K = 0.052;
fofb_type.bpm_sector = 'm1m2c2c3';
fofb_type.corr_sector = 'm1m2c2c3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = [1 80 81 160];
[bpm_idx, corr_idx] = fofb_idx(fofb_type);
%Normalizing models and remove 'NaN' values
fprintf("Normalizing sysid models...\n");
tic
ncorr=size(M,2);
for i=1:ncorr
    if (i==1)||(i==80)||(i==81)||(i==160) %excluded correctors: 1, 80, 81 and 160
        sys{i} = tf(0,1,'Ts',Ts);
    else
        sys{i} = idtf(sys{i}/dcgain(idtf(sys{i})),'Ts',Ts);
    end
end
fprintf("Elapsed time: %f s\n", toc);
fprintf("Building model...\n");
tic
load rf_column.mat
%RF parameters
%https://wiki-sirius.lnls.br/mediawiki/index.php/Machine:RF_System
sync_freq = 2*pi*2.690e3; %Hz -> rad/s
rf_freq = 2*pi*499.664e6;
alpha_c = 1.6e-4;
eta = eta(bpm_idx,1);
M = M(bpm_idx, corr_idx);
Mrf = [M eta];
%Correction matrix with RF correction line
%https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.25.082801
Mcrf = [pinv(M); ((rf_freq*alpha_c)/dot(eta,eta).*eta)'];
sys = sys(corr_idx);
%sys{end+1} = tf(1,1,Ts); %Rf addition
[P, G, C] = ofbmdl_rf(M, eta, Mcrf, K, sys);