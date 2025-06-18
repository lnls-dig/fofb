clc; clear;
addpath 'machine/SIRIUS/'

respmat_fpath = '/ibira/lnls/labs/gie/MachineStudies/FOFBSysId/MATLAB_objs/respmat-no-rf-line-24-05-20.mat';
sysid_res_fpath = '/ibira/lnls/labs/gie/MachineStudies/FOFBSysId/MATLAB_objs/sysid_res-662-23-10-02.mat';

M = load(respmat_fpath).mat_d';
A = load(sysid_res_fpath).sys;

n_corr = size(A,1);
excluded_corr = [1 80 81 160];

Ts = A{2}.Ts;
fs = 1/Ts;

tic
for i=1:n_corr
    if ismember(i,excluded_corr)
        A{i} = tf(0,1,'Ts',Ts);
    else
        A{i} = idtf(A{i}/dcgain(idtf(A{i})),'Ts',Ts);
    end
end

opts = bodeoptions;
opts.FreqUnits = 'Hz';
%opts.XLim = {[0,1e4]};

%figure();
%bode(A{:},opts);
%grid on;

%% Equalization filters

% Picked by inspection
ref_corr_idx = 3;

order = 4;
band = 1e4; % Hz
n_freq_samples = 25;
w = linspace(0,pi,n_freq_samples);

N_butt = 3;
% Low-pass filter as the weighting array
[b_butt,a_butt] = butter(N_butt,band/(fs/2));
h_butt = freqz(b_butt,a_butt,w);
Wt = abs(h_butt);
Wt(w < band/(fs/2)) = 2*Wt(w < band/(fs/2));

F = cell(n_corr,1);
for i=1:n_corr
  if ismember(i,excluded_corr)
    F{i} = tf(0,1,'Ts',Ts);
  elseif i == ref_corr_idx
    F{i} = tf(1,1,'Ts',Ts);
  else
    sys = minreal(absorbDelay(A{ref_corr_idx})/absorbDelay(A{i}));
    %sys = absorbDelay(A{ref_corr_idx})/absorbDelay(A{i});

    %% Yulewalk: disgusting
    %h = freqz(sys.Numerator{1},sys.Denominator{1},f,fs);
    %h(not_interest) = 0;
    %[b,a] = yulewalk(order,f/(fs/2),abs(h));

    h = freqz(sys.Numerator{1},sys.Denominator{1},w);
    %h = h.*h_butt;
    %h = h.*abs(h_butt);
    [b_sys,a_sys] = invfreqz(h,w,order,order,Wt,1000);
    F{i} = tf(b_sys,a_sys,Ts);
    F{i} = F{i}/dcgain(F{i});

    %balredopts = balredOptions;
    %balredopts.StateProjection = 'Truncate';
    %balredopts.FreqIntervals = [0,2*pi*band];
    %F{i} = balred(sys,order,balredopts);

    assert(isstable(F{i}), 'Filter for corrector %s is not stable :c', ...
      int2str(i));
    assert(isstable(A{i}*F{i}));
  end
end

save('F_eq.mat','F');

%figure();
%bode(F{:},opts);
%grid on;
%
%AF = cell(n_corr,1);
%for i=1:n_corr
%  AF{i} = A{i}*F{i};
%end
%
%figure();
%bode(AF{:},opts);
%grid on;

%%

fofb_type.bpm_sector = 'M1M2C2C3';
fofb_type.corr_sector = 'M1M2C2C3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = excluded_corr;
[bpm_idx,corr_idx] = fofb_idx(fofb_type);

M = M(bpm_idx,corr_idx);
Mc = pinv(M);
K = (1/Ts)*[0.12*ones(size(Mc,1)/2,1); 0.166*ones(size(Mc,1)/2,1)];
A = A(corr_idx);
F = F(corr_idx);

figure();

sigmaopts = sigmaoptions;
sigmaopts.FreqUnits = 'Hz';

[P,G,C] = ofbmdl(M,Mc,K,A,[],[],[],[]);
assert(isstable(P))

sigmaplot(P('yd','d'),sigmaopts);
hold on;

[P,G,C] = ofbmdl(M,Mc,K,A,F,[],[],[]);
isstable(P)

sigmaplot(P('yd','d'),sigmaopts);

[P,G,C] = ofbmdl(M,Mc,K,A{ref_corr_idx},[],[],[],[]);
assert(isstable(P))

sigmaplot(P('yd','d'),sigmaopts);
grid on;

legend('Unmatched',strcat('Match attempt (',int2str(order),')'),'Matched');

fprintf('Elapsed time: %.2f s\n',toc);
