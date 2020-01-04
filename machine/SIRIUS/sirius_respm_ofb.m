if ~exist('THERING', 'var')
    sirius;
    global THERING;
end

% Find relevant accelerator components
famdata = sirius_si_family_data(THERING);
bpm = famdata.BPM.ATIndex;
hcm = famdata.CH.ATIndex;
vcm = famdata.CV.ATIndex;
fcm = famdata.FC.ATIndex;

[hcm_fcm, hcm_fcm_sort] = sort([hcm; fcm]);
hfcm_idx = find(hcm_fcm_sort > length(hcm));

[vcm_fcm, vcm_fcm_sort] = sort([vcm; fcm]);
vfcm_idx = find(vcm_fcm_sort > length(vcm));

bpm_fofb = sort([bpm(1:8:end); bpm(4:8:end); bpm(5:8:end); bpm(8:8:end)]);
bpm_sofb = setdiff(bpm, bpm_fofb);

hcm_fofb = fcm;
vcm_fofb = fcm;
hcm_sofb = hcm;
vcm_sofb = vcm;

% Retrieve accelerator parameters: Beta, Mu, Eta, transversal tunes,
% momentum compaction factor (alpha_c) and length (L)
nelm = length(THERING);

[TD, tune, chrom] = twissring(THERING, 0, 1:nelm+1, 'chrom', 1e-8);
beta = reshape([TD.beta]', 2, nelm+1)';
mu = reshape([TD.mu]', 2, nelm+1)';
eta = [TD.Dispersion]';

alpha_c = mcf(THERING);
L = findspos(THERING, nelm+1);

cav_index = findcells(THERING, 'Frequency');
frev = THERING{cav_index}.Frequency/THERING{cav_index}.HarmNumber;
tau = 20e-3*frev;

% Calculate theoretical orbit response matrices
Mh_sofb = respmtheor(beta(:,1), mu(:,1), eta(:,1), alpha_c, L, tune(1), bpm_sofb, hcm_sofb);
Mv_sofb = respmtheor(beta(:,2), mu(:,2), eta(:,3), alpha_c, L, tune(2), bpm_sofb, vcm_sofb);
Mh_fofb = respmtheor(beta(:,1), mu(:,1), eta(:,1), alpha_c, L, tune(1), bpm_fofb, hcm_fofb, 'dynamic', tau);
Mv_fofb = respmtheor(beta(:,2), mu(:,2), eta(:,3), alpha_c, L, tune(2), bpm_fofb, vcm_fofb, 'dynamic', tau);