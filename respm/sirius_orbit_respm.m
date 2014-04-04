function [M,Mrf,Mid,Mdq,Mdbo,Mdbc,Mdbi] = sirius_orbit_respm(THERING, deltakick, deltafreq, displacement)
%SIRIUS_ORBIT_RESPM Orbit response matrices and disturbance matrices.
%   [M,Mrf,Mid,Mdq,Mdbo,Mdbc,Mdbi] = sirius_orbit_respm(THERING, 
%   deltakick, deltafreq, deltakick, displacement) calculates Sirius's 
%   orbit response matrices in 4-D transverse phase space:
%
%   M       corrector magnet kicks inputs [rad]
%   Mrf     RF frequency steps inputs [Hz]
%   Mid     insertion device residual dipolar field disturbance [rad]
%   Mdq     quadrupoles displacements [m]
%   Mdbo    BO bending magnet displacements [m]
%   Mdbc    BC bending magnet displacements [m]
%   Mdbi    BI bending magnet displacements [m]
%
%   Inputs:
%       THERING         AT accelerator model under analysis
%       deltakick       amplitude of dipolar field kicks (rad)
%       deltafreq       amplitude of RF frequency step (Hz)
%       displacement    amplitude of magnets displacement (m)
%
%   Outputs:
%       (Matrices M, Mrf, Mid, Mdq, Mdbo, Mdbc, Mdbi)
%
%   All response matrices are 3-D matrices where each dimension has the
%   following meaning:
%       Dim 1: index of beam orbit value in a given ring position
%       Dim 2: index of input or disturbance affecting the beam orbit
%       Dim 3: selection of 4-D transverse phase space variable where: 
%              1 = horizontal beam position [m]
%              2 = horizontal beam divergence [rad]
%              3 = vertical beam position [m]
%              4 = vertical beam divergence [rad]
%
%   *** FIXME: must include explanation about chosen orbit points; return
%   orbit point indexes, input and disturbance indexes ***

fprintf('\n   -------------------------------\n    Starting "sirius_orbit_respm"\n   -------------------------------\n');

% Dipoles segmentation numbers
nsegs_bo = 32;
nsegs_bc = 44;
nsegs_bi = 16;

%% Markers
bpm = findcells(THERING, 'FamName', 'bpm');

hcm = findcells(THERING, 'FamName', 'cm');
vcm = hcm;

ms  = findcells(THERING, 'FamName', 'ms');  % short straight section
mm  = findcells(THERING, 'FamName', 'mm');  % medium straight section
ml  = findcells(THERING, 'FamName', 'ml');  % long straight section
mc  = findcells(THERING, 'FamName', 'mc');  % bending magnet center

quads = findcells(THERING, 'K');            % quads and dipoles
bends = findcells(THERING, 'BendingAngle');
quads = setdiff(quads, bends);

rf = findcells(THERING, 'FamName', 'cav');

bo = findcells(THERING, 'FamName', 'bo');
bc = findcells(THERING, 'FamName', 'bc');
bi = findcells(THERING, 'FamName', 'bi');

bo = reshape(bo, nsegs_bo, []);
bc = reshape(bc, nsegs_bc, []);
bi = reshape(bi, nsegs_bi, []);

n_bo_points = size(bo,2);
n_bc_points = size(bc,2);
n_bi_points = size(bi,2);

id_points = sort([ms mm ml]);
light_source_points =  sort([mc id_points]);
orbit_points = sort([bpm light_source_points]);

n_id_points = length(id_points);
n_light_source_points = length(light_source_points);
n_orbit_points = length(orbit_points);


%% Calculating response matrix (orbit vs. corrector kicks)
fprintf('   Calculating response matrix (orbit vs. corrector magnet kicks)\n'); tic;
M = zeros(n_orbit_points, length(hcm)+length(vcm), 4);
for i=1:length(hcm)
    original_kick = THERING{hcm(i)}.KickAngle(1);

    THERING{hcm(i)}.KickAngle(1) = original_kick - (deltakick/2);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{hcm(i)}.KickAngle(1) = original_kick + (deltakick/2);
    orbit2 = findorbit4(THERING, 0, orbit_points);

    THERING{hcm(i)}.KickAngle(1) = original_kick;
    
    M(:,i,:)  = (orbit2 - orbit1)'/deltakick;
end
for j=1:length(vcm)
    original_kick = THERING{vcm(j)}.KickAngle(2);
    
    THERING{vcm(j)}.KickAngle(2) = original_kick - (deltakick/2);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{vcm(j)}.KickAngle(2) = original_kick + (deltakick/2);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    
    THERING{vcm(j)}.KickAngle(2) = original_kick;
    
    M(:,i+j,:) = (orbit2 - orbit1)'/deltakick;
end
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);


%% Calculating response matrix (orbit vs. RF frequency)
fprintf('   Calculating response matrix (orbit vs. RF frequency)\n'); tic;
original_cavity_status = getcavity;
setcavity('on');
setradiation('on');

Mrf = zeros(n_orbit_points, 1, 4);
orbit1 = findorbit6(THERING, orbit_points);
THERING{rf}.Frequency = THERING{rf}.Frequency + deltafreq;
orbit2 = findorbit6(THERING, orbit_points);
THERING{rf}.Frequency = THERING{rf}.Frequency - deltafreq;
Mrf(:,1,:)  = (orbit2(1:4,:) - orbit1(1:4,:))'/deltafreq;

setradiation('off');
setcavity(original_cavity_status);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

%% Calculating response matrix (orbit vs. ID center dipole disturbance)
fprintf('   Calculating response matrix (orbit vs. ID center dipole disturbance)\n'); tic;
Mid = zeros(n_orbit_points, 2*n_id_points, 4);
for i=1:length(n_id_points)
    THERING{id_points(i)}.PassMethod = 'CorrectorPass';
    
    THERING{id_points(i)}.KickAngle(1) = deltakick/2;
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{id_points(i)}.KickAngle(1) = -deltakick/2;
    orbit2 = findorbit4(THERING, 0, orbit_points);
    
    THERING{id_points(i)}.KickAngle(1) = 0;
    
    Mid(:,i,:)  = (orbit2 - orbit1)'/deltakick;
    
    THERING{id_points(i)}.KickAngle(2) = deltakick/2;
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{id_points(i)}.KickAngle(2) = -deltakick/2;
    orbit2 = findorbit4(THERING, 0, orbit_points);
    
    THERING{id_points(i)}.KickAngle(2) = 0;
    
    Mid(:,i,:)  = (orbit2 - orbit1)'/deltakick;
    
    THERING{id_points(i)}.PassMethod = 'IdentityPass';
end
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

%% Calculating response matrix (orbit vs. dipole and quadrupole displacements) 
fprintf('   Calculating response matrix (orbit vs. quadrupoles)\n'); tic;
Mdq = displacement_respm(THERING, orbit_points, quads, displacement);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

fprintf('   Calculating response matrix (orbit vs. BO dipoles)\n'); tic;
Mdbo = displacement_respm(THERING, orbit_points, bo, displacement);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

fprintf('   Calculating response matrix (orbit vs. BC dipoles)\n'); tic;
Mdbc = displacement_respm(THERING, orbit_points, bc, displacement);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

fprintf('   Calculating response matrix (orbit vs. BI dipoles)\n'); tic;
Mdbi = displacement_respm(THERING, orbit_points, bi, displacement);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);