% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
%Input data from experiment
%Change ofbmdl.m to reflect rfline status
chosen_bpm = 19;
y_data = h5read('rf-phase-fofb-on-rfline.h5','/data/orbx');
sampling_frequency = h5read('rf-phase-fofb-off.h5','/data/sampling_frequency');
Ts=1/sampling_frequency;
fir_moving_average = dfilt.dffir(ones(1,4)/4);
%197425
y = double(y_data(chosen_bpm,479075:479345));
y = detrend(y,0);
y = filter(fir_moving_average, y);
%y=y*0.65e3;
t = linspace(0,length(y)*Ts,length(y));
hold on;
plot(t,y,'red')

%Creates model
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
fprintf("Elapsed time: %f s\n", toc);
fprintf("Building model...\n");
tic
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
load rf_column.mat
[P, G, C] = ofbmdl_limited_C(M, Mc, K, sys, [],[],[],[],eta(bpm_idx),true);
fprintf("Elapsed time: %f s\n", toc);

%RF step with a closed FOFB
fprintf("RF step - Closed FOFB...\n");
hold on;
for i=19:10:80
    fprintf("Output %d\n",i);
    SYS = P('yd','rf');
    step(SYS(i,1))
end
%figure;
%SYS = P('yd','u(157)');
%step(SYS(chosen_bpm,1))
%hold on;
%for i=1:160
%    fprintf("Output %d\n",i);
%    SYS = P('yd','u(157)');
%    step(SYS(i,1))
%end