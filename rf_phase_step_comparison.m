% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc; clear;
addpath 'machine/SIRIUS/'

respmat_fpath = '/ibira/lnls/labs/gie/MachineStudies/FOFBSysId/MATLAB_objs/respmat-no-rf-line-24-05-20.mat';
sysid_res_fpath = '/ibira/lnls/labs/gie/MachineStudies/FOFBSysId/MATLAB_objs/sysid_res-24-08-13.mat';
rf_column_fpath = '/ibira/lnls/labs/gie/MachineStudies/FOFBSysId/MATLAB_objs/rf_column-24-05-20.mat';

%% FOFB model params
M = load(respmat_fpath).mat_d';
eta = load(rf_column_fpath).eta;
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
eta = eta(bpm_idx)*1.75;
K = (1/Ts)*[0.12*ones(size(M,2)/2,1); 0.166*ones(size(M,2)/2,1)];
A = A(corr_idx);

mov_avg_filt = dfilt.dffir(ones(1,4)/4);
bpm = 19;

%% With RF column

% Creates model
Mc = pinv([M eta]);
Mc(157,:) = []; % discard RF line

[P,G,C] = ofbmdl(M,Mc,K,A,[],[],[],[],eta);

% Data from machine study
orbx = h5read('/ibira/lnls/labs/swc/MachineStudies/22-07-2024/rf-phase-fofb-on-rfline.h5','/data/orbx');
y = filter(mov_avg_filt,double(orbx(bpm,:)));
y = y(140:339);
y = detrend(y,0);
y = y*1e3; % adjusting to nm
t = linspace(0,(length(y)-1)*Ts,length(y));

% TODO: phase adjustment
y = -y;

figure;
hold on;
plot(t,y,'Color','#D95319');
sys = P('yd','ph');
step(sys(bpm,157),t);
title('With RF column');
ylabel('Position displacement [nm]')
legend({'Data','Model'});
grid on;

%% Without RF column
% Creates model
Mc = pinv(M);
[P,G,C] = ofbmdl(M,Mc,K,A,[],[],[],[],eta);

% Data from machine study
orbx = h5read('/ibira/lnls/labs/swc/MachineStudies/22-07-2024/rf-phase-fofb-on.h5','/data/orbx');
y = filter(mov_avg_filt,double(orbx(bpm,:)));
y = y(5526:5724);
y = detrend(y,0);
y = y*1e3; % adjusting to nm
t = linspace(0,(length(y)-1)*Ts,length(y));

% TODO: phase adjustment
y = -y;

figure;
hold on;
plot(t,y,'Color','#D95319');
sys = P('yd','ph');
step(sys(bpm,157),t);
title('Without RF column');
ylabel('Position displacement [nm]')
legend({'Data','Model'});
grid on;