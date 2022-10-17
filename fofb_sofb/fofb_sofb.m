%%

clear

plane = 'h';
include_rf = true;

nbpm_sec = 8;
nch_sec = 6;
ncv_sec = 8;
nfc_sec = 4;
nsec_total = 20;

nsec = 4;

%%
nbpms = nbpm_sec*nsec;
nchs = nch_sec*nsec;
ncvs = ncv_sec*nsec;
nfcs = nfc_sec*nsec;

nbpms_total = nbpm_sec*nsec_total;
nchs_total = nch_sec*nsec_total;
ncvs_total = ncv_sec*nsec_total;
nfcs_total = nfc_sec*nsec_total;

if strcmpi(plane, 'h')
    bpm_idx = [1:nbpms-1 nbpms_total];
    corr_idx_s = [1:nchs-1 nchs_total];
    corr_idx_f = [1:nfcs-1 nfcs_total];
    sel_bpm_f = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1]) length(bpm_idx)];
    sel_corr_f = [sort([nfc_sec:nfc_sec:nfcs-1 1:nfc_sec:nfcs-1]) length(corr_idx_f)];
elseif strcmpi(plane, 'v')
    bpm_idx = nbpms_total + [1:nbpms-1 nbpms_total];
    corr_idx_s = nchs_total + [1:ncvs-1 ncvs_total];
    corr_idx_f = nfcs_total + [1:nfcs-1 nfcs_total];
    sel_bpm_f = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1]) length(bpm_idx)];
    sel_corr_f = [sort([nfc_sec:nfc_sec:nfcs-1 1:nfc_sec:nfcs-1]) length(corr_idx_f)];
end
sel_posang = sort([1:nbpm_sec/2:nbpms/2 3:nbpm_sec/2:nbpms/2]);
sel_posang = [sel_posang nbpms_total/2+sel_posang];

if include_rf
    corr_idx_s = [corr_idx_s nchs_total+ncvs_total+1];
    corr_idx_f = [corr_idx_f 2*nfcs_total+1];
    
    sel_corr_f = [sel_corr_f nfcs+1];
end

[Ms,Mcs] = fofb_sofb_matrices('sofb_respmat.txt', bpm_idx, corr_idx_s, [], []);
[Mf,Mcf] = fofb_sofb_matrices('fofb_respmat.txt', bpm_idx, corr_idx_f, sel_bpm_f, sel_corr_f);
[posx, posy, angx, angy] = bpmpos2posang_matrix;

if strcmpi(plane, 'h')
    Wz = [posx.bpmx(:, bpm_idx); angx.bpmx(:, bpm_idx)];
elseif strcmpi(plane, 'v')
    Wz = [posy.bpmy(:, bpm_idx); angy.bpmy(:, bpm_idx)];
end
Wz = Wz(sel_posang,:);

%%
% dly [s]
% bw  [Hz]

fofb_sofb_param(1).M = Ms;
fofb_sofb_param(1).Mc = Mcs;
fofb_sofb_param(1).dly = 80e-3;
fofb_sofb_param(1).H = tf(1, [1/2/pi/10 1]);
fofb_sofb_param(1).Ki = 0.5;
fofb_sofb_param(1).char = 's';
fofb_sofb_param(1).bpm_idx = bpm_idx;
fofb_sofb_param(1).corr_idx = corr_idx_s;
fofb_sofb_param(1).sel_bpm = [];
fofb_sofb_param(1).sel_corr = [];

fofb_sofb_param(2).M = Mf;
fofb_sofb_param(2).Mc = Mcf;
fofb_sofb_param(2).dly = 25e-6;
fofb_sofb_param(2).H = tf(1, [1/2/pi/10e3 1]);
fofb_sofb_param(2).Ki = 5000;
fofb_sofb_param(2).char = 'f';
fofb_sofb_param(2).bpm_idx = bpm_idx;
fofb_sofb_param(2).corr_idx = corr_idx_f;
fofb_sofb_param(2).sel_bpm = sel_bpm_f;
fofb_sofb_param(2).sel_corr = sel_corr_f;

acfofb_sofb_param = fofb_sofb_param;
acfofb_sofb_param(2).H = fofb_sofb_param(2).H*tf([1 0],[1 2*pi*1]);

sofb_param = fofb_sofb_param;
sofb_param(2).Ki = 0;

fofb_param = fofb_sofb_param;
fofb_param(1).Ki = 0;
fofb_param(2).sel_bpm = [];
fofb_param(2).sel_corr = [];

param = {fofb_sofb_param acfofb_sofb_param sofb_param fofb_param};
param_names = {'fofb_sofb', 'acfofb_sofb', 'sofb', 'fofb'};

%% Compute transfer functions
for i=1:length(param)
    T{i} = fofb_sofb_tf(param{i}, [], Wz);
    T{i}.name = param_names{i};
    d_us{i} = getIOTransfer(T{i},'d','us');
    d_uf{i} = getIOTransfer(T{i},'d','uf');
    d_ydd{i} = getIOTransfer(T{i},'d','ydd');
    d_z{i} = getIOTransfer(T{i},'d','z');
    
    d_us{i}.name = [param_names{i} ' d->us'];
    d_uf{i}.name = [param_names{i} ' d->uf'];
    d_ydd{i}.name = [param_names{i} ' d->ydd'];

    d_z{i}.name = [param_names{i} ' d->z'];
end


%%
figure
sigma(d_ydd{:}, 2*pi*logspace(-5, 4, 1e3));
legend show
