function [bpm_idx, corr_idx] = fofb_idx(fofb_type)
% fofb_type = 'm1m2' or 'm1m2c2c3'

nbpm_sec = 8;
ncorr_sec = 4;
nsec_total = 20;

%%
nbpms = nbpm_sec*nsec_total;
bpm_idx_h = 1:nbpms;
bpm_idx_v = nbpms + bpm_idx_h;

ncorrs = ncorr_sec*nsec_total;
corr_idx_h = 1:ncorrs;
corr_idx_v = ncorrs + corr_idx_h;

if strcmpi(fofb_type.bpm_sector, 'm1m2')
    sel_bpm = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1]) nbpms];
elseif strcmpi(fofb_type.bpm_sector, 'm1m2c2c3')
    sel_bpm = [sort([nbpm_sec:nbpm_sec:nbpms-1 1:nbpm_sec:nbpms-1 4:nbpm_sec:nbpms-1 5:nbpm_sec:nbpms-1]) nbpms];
end

if strcmpi(fofb_type.corr_sector, 'm1m2')
    sel_corr = [sort([ncorr_sec:ncorr_sec:ncorrs-1 1:ncorr_sec:ncorrs-1]) ncorrs];
elseif strcmpi(fofb_type.corr_sector, 'm1m2c2c3')
    sel_corr = [sort([ncorr_sec:ncorr_sec:ncorrs-1 1:ncorr_sec:ncorrs-1 2:ncorr_sec:ncorrs-1 3:ncorr_sec:ncorrs-1]) ncorrs];
end

bpm_idx = [bpm_idx_h bpm_idx_v];
corr_idx = [corr_idx_h corr_idx_v];

if isfield(fofb_type, 'bpm_remove_idx')
    bpm_idx = bpm_idx(setdiff([sel_bpm nbpms+sel_bpm], fofb_type.bpm_remove_idx));    
end

if isfield(fofb_type, 'corr_remove_idx')
    corr_idx = corr_idx(setdiff([sel_corr ncorrs+sel_corr], fofb_type.corr_remove_idx));
end

