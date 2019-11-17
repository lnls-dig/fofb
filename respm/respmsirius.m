function [Mcorr, Mss, Mdisp, Mrf] = respmsirius(plane)
%RESPMSIRIUS Orbit response matrices and disturbance matrices.
%   M = fofb_orbit_respm_sirius(THERING) calculates Sirius's orbit response
%   matrices in 4-D transverse phase space:
%
%           corrector magnet kicks inputs [rad]
%           RF frequency steps inputs [Hz]
%           insertion device residual dipolar field disturbance [rad]
%           quadrupoles displacements [m]
%           bending magnets displacements [m]
%
%   Inputs:
%       THERING: AT accelerator model of Sirius ring
%
%   Outputs:
%       M: 
%
%   All response matrices are 3-D matrices where each dimension has the
%   following meaning:
%       Dim 1: index of beam orbit value in a given ring position
%       Dim 2: index of input or disturbance affecting the beam orbit
%       Dim 3: 4-D transverse phase space variable: 
%              1 = horizontal beam position [m]
%              2 = horizontal beam angle [rad]
%              3 = vertical beam position [m]
%              4 = vertical beam angle [rad]
%
%   *** FIXME: must include explanation about chosen orbit points; return
%   orbit point indexes, input and disturbance indexes ***

global THERING;
% sirius;
% setoperationalmode(1);

fprintf('\n   -------------------------------\n    Starting "sirius_orbit_respm"\n   -------------------------------\n');

% Dipoles segmentation numbers
nsegs_bc = 36;
nsegs_b1 = 32;
nsegs_b2 = 38;

% Markers

% Insertion devices
mc = findcells(THERING, 'FamName', 'mc')';
mia  = findcells(THERING, 'FamName', 'mia')';
mib  = findcells(THERING, 'FamName', 'mib')';
mip  = findcells(THERING, 'FamName', 'mip')';

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

% Dedicated vertical correctors
FC1 = findcells(THERING, 'FamName', 'FC1')';
FC2 = findcells(THERING, 'FamName', 'FC2')';

% Dipoles
BC = findcells(THERING, 'FamName', 'BC')';
B1 = findcells(THERING, 'FamName', 'B1')';
B2 = findcells(THERING, 'FamName', 'B2')';
BC_edge = findcells(THERING, 'FamName', 'BC_EDGE')';
B1_edge = findcells(THERING, 'FamName', 'B1_edge')';
B2_edge = findcells(THERING, 'FamName', 'B2_edge')';

BC = sort([BC; BC_edge]);
B1 = sort([B1; B1_edge]);
B2 = sort([B2; B2_edge]);

% Takes dipole segmentation into account
BC = reshape(BC, nsegs_bc, []);
B1 = reshape(B1, nsegs_b1, []);
B2 = reshape(B2, nsegs_b2, []);

bpm = findcells(THERING, 'FamName', 'BPM')';
id = sort([mia; mib]);
source =  sort([mc; id]);
hcm = sort([SDA0; SDA1; SDB1; SFB0; SFP0; SDP1; SFP2; SFA2; SFB2_CH]);
vcm = sort([SDA0; SDA1; SDB1; SFB0; SFP0; SDP1; SFP2; SDA3; SDB3; SDP3; SFB2_CV; CV]);
crhv = sort([FC1; FC2]);
quad = sort([QDA; QFA; QDB1; QDB2; QFB; Q1; Q2; Q3; Q4]);
%dipole = sort(BC B1 B2])';
rf = findcells(THERING, 'FamName', 'cav')';

if strcmpi(plane, 'h')
    cm = hcm;
elseif strcmpi(plane, 'v')
    cm = vcm;
end

markers = struct( ...
    'orbit', bpm, ...
    'corr', vcm ...
);
%    'ss', id, ...
%    'disp', quad, ...
%    'rf', rf ...
%);

[Mcorr, Mss, Mdisp, Mrf] = respmorbit(THERING, markers, plane);
