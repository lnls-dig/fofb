% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
%FOFB parameters
load respmat.mat
M = mat_d';
sys=cell(size(M,2), 1);
load sysid_res-24-08-13.mat
Ts=1/48193;
K = 0.052*48193;
fofb_type.bpm_sector = 'm1m2c2c3';
fofb_type.corr_sector = 'm1m2c2c3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = [1 80 81 160];
[bpm_idx, corr_idx] = fofb_idx(fofb_type);
%Normalizing models and removing 'NaN' values
fprintf("Normalizing sysid models...\n");
tic
for i=1:size(M,2)
    if (i==1)||(i==80)||(i==81)||(i==160) %excluded correctors: 1, 80, 81 and 160
        sys{i} = tf(0,1,'Ts',Ts);
    else
        sys{i} = idtf(sys{i}/dcgain(idtf(sys{i})),'Ts',Ts);
    end
end

for i=1:size(M,2)
    if (i==1)||(i==80)||(i==81)||(i==160) %excluded correctors: 1, 80, 81 and 160
        sys{i} = tf(0,1,'Ts',Ts);
    else
        sys{i} = sys{5};
    end
end

fprintf("Elapsed time: %f s\n", toc);
fprintf("Building model...\n");
tic
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
load rf_column.mat
[P, G, C] = ofbmdl(M, Mc, K, sys, [],[],[],[],eta(bpm_idx));
fprintf("Elapsed time: %f s\n", toc);
fofb.P = P;
fofb.G = G;
fofb.C = C;
fofb.S = P('yd','d');
fofb.S_noise = P('yd',{'d','n'});
%Save FOFb model
save('fofb_model_rf','fofb');

%RF step with a closed FOFB
fprintf("RF step - Closed FOFB...\n");
figure;
hold on;
for i=1:160
    fprintf("Output %d\n",i);
    SYS = P('yd','u(157)');
    step(SYS(i,1))
end