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

% Markers
family_data = sirius_si_family_data(THERING);

% Insertion devices
mc = findcells(THERING, 'FamName', 'mc')';
mia  = findcells(THERING, 'FamName', 'mia')';
mib  = findcells(THERING, 'FamName', 'mib')';
mip  = findcells(THERING, 'FamName', 'mip')';

% Quadrupoles
qfa  = findcells(THERING, 'FamName', 'QFA')';
qfb  = findcells(THERING, 'FamName', 'QFB')';
qfp  = findcells(THERING, 'FamName', 'QFP')';
qda = findcells(THERING, 'FamName', 'QDA')';
qdb1 = findcells(THERING, 'FamName', 'QDB1')';
qdb2 = findcells(THERING, 'FamName', 'QDB2')';
qdp1  = findcells(THERING, 'FamName', 'QDP1')';
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
bpm = sort(family_data.BPM.ATIndex);

% Orbit correctors
fhcm = sort(family_data.FCH.ATIndex);
fvcm = sort(family_data.FCV.ATIndex);
hcm = sort(family_data.CH.ATIndex);
vcm = sort(family_data.CV.ATIndex);

bpm = findcells(THERING, 'FamName', 'BPM')';
id = sort([mia; mib; mip]);
source =  sort([mc; id]);
quad = sort([qda; qfa; qdb1; qdb2; qfb; qf1; qf2; qf3; qf4; qfp; qdp1; qdp2]);
dipole = sort([bc b1 b2])';
rf = findcells(THERING, 'FamName', 'SRFCav')';

markers = struct( ...
    'orbit', bpm, ...
    'corr', hcm, ...
    'ss', id, ...
    'disp', quad, ...
    'rf', rf ...
);

[Mcorr, Mss, Mdisp, Mrf] = respmorbit(THERING, markers, plane);