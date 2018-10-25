function r = respmsirius
%RESPMSIRIUS Orbit response matrices and disturbance matrices.
%   r = fofb_orbit_respm_sirius

global THERING;

% Markers
family_data = sirius_si_family_data(THERING);

% Insertion devices
mc  = findcells(THERING, 'FamName', 'mc')';
mia = findcells(THERING, 'FamName', 'mia')';
mib = findcells(THERING, 'FamName', 'mib')';
mip = findcells(THERING, 'FamName', 'mip')';

% Quadrupoles
qfa  = findcells(THERING, 'FamName', 'QFA')';
qfb  = findcells(THERING, 'FamName', 'QFB')';
qfp  = findcells(THERING, 'FamName', 'QFP')';
qda  = findcells(THERING, 'FamName', 'QDA')';
qdb1 = findcells(THERING, 'FamName', 'QDB1')';
qdb2 = findcells(THERING, 'FamName', 'QDB2')';
qdp1 = findcells(THERING, 'FamName', 'QDP1')';
qdp2 = findcells(THERING, 'FamName', 'QDP2')';
qf1  = findcells(THERING, 'FamName', 'Q1')';
qf2  = findcells(THERING, 'FamName', 'Q2')';
qf3  = findcells(THERING, 'FamName', 'Q3')';
qf4  = findcells(THERING, 'FamName', 'Q4')';

% Dipoles
bc = sort([findcells(THERING, 'FamName', 'BC') findcells(THERING, 'FamName', 'BC_EDGE')]);
b1 = sort([findcells(THERING, 'FamName', 'B1') findcells(THERING, 'FamName', 'B1_EDGE')]);
b2 = sort([findcells(THERING, 'FamName', 'B2') findcells(THERING, 'FamName', 'B2_EDGE')]);

% Takes dipole segmentation into account
bc = reshape(bc, [], 20)';
b1 = reshape(b1, [], 20)';
b2 = reshape(b2, [], 20)';

% BPMs
bpm_id = sort([findcells(THERING, 'FamName', 'id_enda') findcells(THERING, 'FamName', 'id_endb') findcells(THERING, 'FamName', 'id_endp')])';
bpm_ring = sort(family_data.BPM.ATIndex);
bpm = sort([bpm_id; bpm_ring]);

% Orbit correctors
hcm_slow = sort(family_data.CH.ATIndex);
vcm_slow = sort(family_data.CV.ATIndex);
hcm_fast = sort([family_data.FC1.ATIndex; family_data.FC2.ATIndex]);
vcm_fast = hcm_fast;

id = sort([mia; mib; mip]);
light_source =  sort([mc; id]);
quad = sort([qda; qfa; qdb1; qdb2; qfb; qf1; qf2; qf3; qf4; qfp; qdp1; qdp2]);
%dipole = sort([bc b1 b2])';

[r.bpm.sofb.Mx, r.bpm.sofb.Mx_] = respmcorr(bpm, hcm_slow, 'x');
[r.bpm.sofb.My, r.bpm.sofb.My_] = respmcorr(bpm, vcm_slow, 'y');
[r.bpm.fofb.Mx, r.bpm.fofb.Mx_] = respmcorr(bpm, hcm_fast, 'x');
[r.bpm.fofb.My, r.bpm.fofb.My_] = respmcorr(bpm, vcm_fast, 'y');
[r.bpm.dist.id.Mx, r.bpm.dist.id.Mx_] = respmdisp(bpm, id, 'x');
[r.bpm.dist.id.My, r.bpm.dist.id.My_] = respmdisp(bpm, id, 'y');
[r.bpm.dist.misalign.quad.Mx, r.bpm.dist.misalign.quad.Mx_] = respmdisp(bpm, quad , 'x');
[r.bpm.dist.misalign.quad.My, r.bpm.dist.misalign.quad.My_] = respmdisp(bpm, quad , 'y');
[r.bpm.dist.rf.M, r.bpm.dist.rf.M_] = respmrf (bpm);

[r.light_source.sofb.Mx, r.light_source.sofb.Mx_] = respmcorr(light_source, hcm_slow, 'x');
[r.light_source.sofb.My, r.light_source.sofb.My_] = respmcorr(light_source, vcm_slow, 'y');
[r.light_source.fofb.Mx, r.light_source.fofb.Mx_] = respmcorr(light_source, hcm_fast, 'x');
[r.light_source.fofb.My, r.light_source.fofb.My_] = respmcorr(light_source, vcm_fast, 'y');
[r.light_source.dist.id.Mx, r.light_source.dist.id.Mx_] = respmdisp(light_source, id, 'x');
[r.light_source.dist.id.My, r.light_source.dist.id.My_] = respmdisp(light_source, id, 'y');
[r.light_source.dist.misalign.quad.Mx, r.light_source.dist.misalign.quad.Mx_] = respmdisp(light_source, quad , 'x');
[r.light_source.dist.misalign.quad.My, r.light_source.dist.misalign.quad.My_] = respmdisp(light_source, quad , 'y');
[r.light_source.dist.rf.M, r.light_source.dist.rf.M_] = respmrf (light_source);

% fprintf('\n   -------------------------------\n    Starting "sirius_orbit_respm"\n   -------------------------------\n');
% % Response matrix: orbit vs. corrector kicks
%         fprintf('   Response matrix: orbit vs. corrector magnet kicks...\n'); tic;
%         % Response matrix: orbit vs. RF frequency
%         fprintf('   Response matrix: orbit vs. RF frequency...\n'); tic;
%         % Response matrix: orbit vs. magnet's displacements
%         fprintf('   Response matrix: orbit vs. magnets'' displacements...\n'); tic;
%         % Response matrix: orbit vs. ID disturbance
%         fprintf('   Response matrix: orbit vs. ID disturbance...\n'); tic;
% fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);