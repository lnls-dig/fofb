function [Mcorr, Mss, Mdisp, Mrf] = respmsirius(plane)
%RESPMSIRIUS Orbit response matrices and disturbance matrices.
%   r = fofb_orbit_respm_sirius

global THERING;

% Markers

% Insertion devices
mc  = findcells(THERING, 'FamName', 'mc')';
mia = findcells(THERING, 'FamName', 'mia')';
mib = findcells(THERING, 'FamName', 'mib')';
mip = findcells(THERING, 'FamName', 'mip')';

% Quadrupoles
QDA  = findcells(THERING, 'FamName', 'QDA')';
QFA  = findcells(THERING, 'FamName', 'QFA')';
QFB  = findcells(THERING, 'FamName', 'QFB')';
Q1   = findcells(THERING, 'FamName', 'Q1')';
Q2   = findcells(THERING, 'FamName', 'Q2')';
Q3   = findcells(THERING, 'FamName', 'Q3')';
Q4   = findcells(THERING, 'FamName', 'Q4')';
QDB1 = findcells(THERING, 'FamName', 'QDB1')';
QDB2 = findcells(THERING, 'FamName', 'QDB2')';

% Sextupoles with correctors CH and CV
SDA0 = findcells(THERING, 'FamName', 'SDA0')';
SDA1 = findcells(THERING, 'FamName', 'SDA1')';
SDB1 = findcells(THERING, 'FamName', 'SDB1')';
SFB0 = findcells(THERING, 'FamName', 'SFB0')';
SFP0 = findcells(THERING, 'FamName', 'SFP0')';
SDP1 = findcells(THERING, 'FamName', 'SDP1')';
SFP2 = findcells(THERING, 'FamName', 'SFP2')';

% Sextupoles with correctors CH
SFA2 = findcells(THERING, 'FamName', 'SFA2')';

% Sextupoles with correctors CV
SDA3 = findcells(THERING, 'FamName', 'SDA3')';
SDB3 = findcells(THERING, 'FamName', 'SDB3')';
SDP3 = findcells(THERING, 'FamName', 'SDP3')';

% Sextupoles with CH or CV depending on sector
SFB2 = findcells(THERING, 'FamName', 'SFB2')';
SFB2_CH = SFB2;
SFB2_CV = SFB2(1:2:end);

% Dedicated vertical correctors
CV = findcells(THERING, 'FamName', 'CV')';

% Dedicated fast correctors
FC1 = findcells(THERING, 'FamName', 'FC1')';
FC2 = findcells(THERING, 'FamName', 'FC2')';

% Dipoles
bc = sort([findcells(THERING, 'FamName', 'BC') findcells(THERING, 'FamName', 'BC_EDGE')]);
b1 = sort([findcells(THERING, 'FamName', 'B1') findcells(THERING, 'FamName', 'B1_EDGE')]);
b2 = sort([findcells(THERING, 'FamName', 'B2') findcells(THERING, 'FamName', 'B2_EDGE')]);

% Takes dipole segmentation into account
bc = reshape(bc, [], 20)';
b1 = reshape(b1, [], 20)';
b2 = reshape(b2, [], 20)';


% Orbit correctors

bpm = findcells(THERING, 'FamName', 'BPM')';
id = sort([mia; mib; mip]);
source =  sort([mc; id]);
hcm = sort([SDA0; SDA1; SDB1; SFB0; SFP0; SDP1; SFP2; SFA2; SFB2_CH]);
vcm = sort([SDA0; SDA1; SDB1; SFB0; SFP0; SDP1; SFP2; SDA3; SDB3; SDP3; SFB2_CV; CV]);
crhv = sort([FC1; FC2]);
quad = sort([QDA; QFA; QDB1; QDB2; QFB; Q1; Q2; Q3; Q4]);
%dipole = sort([bc b1 b2])';
rf = findcells(THERING, 'FamName', 'cav')';

markers = struct( ...
    'orbit', bpm, ...
    'corr', hcm, ...
    'ss', id, ...
    'disp', quad ...
);

[Mcorr, Mss, Mdisp, Mrf] = respmorbit(THERING, markers, plane);
