%FOFB model with RF transfer function (rf phase noise -> BPM) implemented
%as disturbance

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
%Normalizing models and removing 'NaN' values
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
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
%RF parameters
%https://wiki-sirius.lnls.br/mediawiki/index.php/Machine:RF_System
load rf_column.mat
eta = eta(bpm_idx,1);
sync_freq = 2*pi*2.690e3; %Hz -> rad/s
rf_freq = 2*pi*499.664e6;
alpha_c = 1.6e-4;
%Rf phase to beam transfer function
%https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.25.082801
H = tf([sync_freq^2/rf_freq/alpha_c 0],[1 2*alpha_c sync_freq^2]);
H = c2d(H, Ts)*1e9; %adjusting for nm 
H = H.*diag(eta);
%Rf freq. to beam position as a simple 2nd order system
res_freq = 2*pi*2.69e3;   %ressonance at synchrotron frequency
csi = 0.0000001;
wn = res_freq/sqrt(1-2*csi^2);
H2 = tf(wn^2,[1 2*csi*wn wn^2]);
bode(H(1,1),H2(1,1));
H2 = c2d(H2, Ts)*1e9; %adjusting for nm 
H2 = H2.*diag(eta);
[P, G, C] = ofbmdl(M, Mc, K, sys, tf(1, 1, Ts), H);
fprintf("Elapsed time: %f s\n", toc);