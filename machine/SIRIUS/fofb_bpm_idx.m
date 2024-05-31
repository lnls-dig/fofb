function bpm_idx = fofb_bpm_idx(fofb_bpms)
% fofb_type = 'm1m2' or 'm1m2c2c3'

nbpm_sec = 8;
nsec_total = 20;

%%
nbpms = nbpm_sec*nsec_total;

bpm_idx_h = 1:nbpms;
bpm_idx_v = nbpms + bpm_idx_h;

if strcmpi(fofb_bpms, 'm1m2')
    sel_bpm = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1]) nbpms];
elseif strcmpi(fofb_bpms, 'm1m2c2c3')
    sel_bpm = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1 4:nbpm_sec:nbpms-1 5:nbpm_sec:nbpms-1]) nbpms];
end

bpm_idx = [bpm_idx_h(sel_bpm) bpm_idx_v(sel_bpm)];