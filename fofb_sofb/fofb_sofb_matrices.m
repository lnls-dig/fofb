function [M,Mc] = fofb_sofb_matrices(filename, bpm_idx, corr_idx, sel_bpm, sel_corr)

% Response matrix
M = dlmread(filename,' ',4,0);
M = M(bpm_idx, corr_idx);
if ~isempty(sel_corr)
    M = M(:, sel_corr);
end

% Correction matrix (pseudo-inverse)
if isempty(sel_bpm)
    Mc = pinv(M);
else
    M_ = M;
    M_(setdiff(1:length(bpm_idx), sel_bpm), :) = 0;
    Mc = pinv(M_);
end

end