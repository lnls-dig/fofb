% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
%Compare FOFB actuation with and without RF column
%Common imported data
load respmat.mat
M = mat_d';
sys=cell(size(M,2), 1);
load sysid_res-24-08-13.mat
Ts=1/48193;
K = 0.13*48193;
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
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
load rf_column.mat

%Without RF column
%Input data from experiment
figure;
chosen_bpm = 19;
y_data = h5read('rf-phase-fofb-on.h5','/data/orbx');
sampling_frequency = h5read('rf-phase-fofb-on.h5','/data/sampling_frequency');
Ts=1/sampling_frequency;
fir_moving_average = dfilt.dffir(ones(1,4)/4);
y = double(y_data(chosen_bpm,197425:197540));
y = detrend(y,0);
y = filter(fir_moving_average, y);
%y=y/(2*pi);
t = linspace(0,length(y)*Ts,length(y));
hold on;
plot(t,y,'red')

%Creates model
fprintf("Building model...\n");
tic
[P, ~, ~] = ofbmdl_limited_C(M,[],K,sys,[],[],[],[],eta(bpm_idx),[]);
fprintf("Elapsed time: %f s\n", toc);

%RF step with a closed FOFB
%Input data from experiment
fprintf("RF step - Closed FOFB...\n");
hold on;
for i=chosen_bpm:80:80
    fprintf("Output %d\n",i);
    SYS = P('yd','rf');
    step(SYS(i,1))
end
SYS = P('yd','rf');
[y_fofb_on,~] = step(SYS(chosen_bpm,1));

%figure;
%hold on;
%for i=1:156
%    fprintf("Control Effort %d\n",i);
%    SYS = getIOTransfer(P,"rf","u");
%    step(SYS(i,1))
%    title('Control Effort - Without RF column')
%end

%With RF column
figure;
y_data = h5read('rf-phase-fofb-on-rfline.h5','/data/orbx');
sampling_frequency = h5read('rf-phase-fofb-on-rfline.h5','/data/sampling_frequency');
Ts=1/sampling_frequency;
fir_moving_average = dfilt.dffir(ones(1,4)/4);
y = double(y_data(chosen_bpm,479075:479345));
y = detrend(y,0);
y = filter(fir_moving_average, y);
%y=y/(2*pi);
t = linspace(0,length(y)*Ts,length(y));
hold on;
plot(t,y,'red')

%Creates model
fprintf("Building model...\n");
tic
[P, G, C] = ofbmdl_limited_C(M,[],K,sys,[],[],[],[],eta(bpm_idx),true);
fprintf("Elapsed time: %f s\n", toc);

%RF step with a closed FOFB
fprintf("RF step - Closed FOFB...\n");
hold on;
for i=chosen_bpm:80:80
    fprintf("Output %d\n",i);
    SYS = P('yd','rf');
    step(SYS(i,1))
end
%figure;
%hold on;
%for i=1:156
%    fprintf("Control Effort %d\n",i);
%    SYS = getIOTransfer(P,"rf","u");
%    step(SYS(i,1))
%    title('Control Effort - With RF column')
%end
SYS = P('yd','rf');
[y_fofb_on_rf,~] = step(SYS(chosen_bpm,1));

%Plot two figure in the same plot
figure;
min_v = min(length(y_fofb_on_rf),length(y_fofb_on));
hold on;
plot(t(1:min_v),y_fofb_on_rf(1:min_v),'blue')
plot(t(1:min_v),y_fofb_on(1:min_v),'red')
legend('With RF column', 'Without RF column');
grid on;