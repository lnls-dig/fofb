% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

%Gain  variation
%FOFB parameters
load respmat.mat
M = mat_d';
sys=cell(size(M,2), 1);
load sysid_res-24-08-13.mat
Ts=1/48193;
K = 0.120*48193;
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
fprintf("Elapsed time: %f s\n", toc);
fprintf("Building model...\n");
tic
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
[P, G, C] = ofbmdl(M, Mc, K, sys);
fprintf("Elapsed time: %f s\n", toc);
%Considering correctos discrepancies
stable = isstable(P);
it=0;
while stable
    it=it+1;
    fprintf("Iteration= %d\n", it);
    fprintf("K= %f\n", K/48193);
    K=K*1.1;
    [P, G, C] = ofbmdl(M, Mc, K, sys);
    stable = isstable(P);
end
%Equal correctors
%If only a SISO restriction is considered 
%bpm_idx = 8;
%corr_idx = 4;