function [P, G, C] = ofbmdl(M, Mc, K, A, F, Wd, Wz, Wn, eta)
% OFBMDL Orbit Feedback Model.
% 
% Construct orbit feedback state-space model.
%
% [P, G, C] = OFBMDL(M, Mc, K, A, F, Wd, Wz, Wn)
%
% INPUTS:
%   M:    Orbit response matrix. Dimensions NBPM x NCORR, where NBPM is
%         the number of BPMs and NCORR is the number of orbit correctors.
%   Mc:   Orbit correction matrix. Dimensions NCORR x NBPM.
%   K:    Controller gain (scalar or array with NCORR elements).
%   A:    Actuator responses (dynamic system object (tf, zpk, ss, etc.) or
%         cell array of dynamic system objects with NCORR elements)
%         encompassing power supply, magnet, vacuum chamber, BPM and
%         network delay dynamics. If only one system is provided all are
%         assumed to have identical responses. 
%   F:    (Optional input) Shaping filters per orbit corrector (dynamic
%         system object (tf, zpk, ss, etc.) or cell array of dynamic system
%         objects with NCORR elements). If only one filter is provided the
%         shaping filters of all actuators are assumed to be identical.
%   Wd:   (Optional input) Disturbance weighting matrix, used for loop
%         optimzation purposes.
%   Wz:   (Optional input) Performance weighting matrix, used for loop
%         optimzation purposes.
%   Wn:   (Optional input) Noise weighting matrix, used for loop
%         optimzation purposes.
%   eta:  (Optional input) RF column in the orbit response matrix. Dimensions NBPMx1.
%
% OUTPUTS:
%   P:    Generalized plant of the closed-loop orbit feedback system.
%   G:    Combined responses of actuators and beam (only static beam
%         response).
%   C:    Combined responses of the controller and shaping filters.

nbpms = size(M,1);
ncorr = size(M,2);
rf_interaction=true;

if nargin < 5 || isempty(F)
    F = tf(1, 1, -1);
end

if nargin < 6 || isempty(Wd)
    Wd = eye(nbpms);
end

if nargin < 7 || isempty(Wz)
    Wz = eye(nbpms);
end

if nargin < 8 || isempty(Wn)
    Wn = eye(nbpms);
end

if nargin < 9 || isempty(eta)
    rf_interaction=false;
end

% Actuator transfer function (includes power supply, magnet, vacuum
% chamber, BPM, beam and network delay dynamics)
[A, Ts] = sys_array2matrix(A, ncorr);

% Shaping filters
[F, Ts_F] = sys_array2matrix(F, ncorr);

if Ts_F ~= -1 && Ts ~= Ts_F
    error('Sample time ''Ts'' of actuator transfer functions and shaping filters must be equal.');
end

%Adds optional RF interaction
if rf_interaction
    %RF parameters obtained from fit_rf.m and wiki-sirius
    %https://wiki-sirius.lnls.br/mediawiki/index.php/Machine:RF_System
    sync_freq = 2.1515e3; %Hz (Fit)
    alpha_s = 1.2547e3; %(Fit)
    rf_freq = 499.664e6; %Hz (Theoretical)
    alpha_c = 1.6e-4; %(Theoretical)
    %RF phase to beam transfer function
    %https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.25.082801
    H = tf([sync_freq^2/rf_freq/alpha_c 0],[1 2*alpha_s sync_freq^2]);
    H = c2d(H, Ts);%*1e9; %adjusting for nm
    H = H.*eta;
    % New M and altered Mc
    Mrf = [M H];
    Mc = pinv([M zeros(length(eta),1)]);
    Arf = [A; ones(ncorr,1)'];
    % Altered shaping filters
    F(ncorr+1,ncorr+1)=tf(1, 1, Ts);
    % Plant transfer function
    G = ss(Mrf)*ss(Arf);
    % Controller transfer function
    C = ss(F)*diag(K)*Mc*tf([Ts 0], [1 -1], Ts);
else
    % Plant transfer function
    G = ss(M)*ss(A);
    % Controller transfer function
    C = ss(F)*diag(K)*Mc*tf([Ts 0], [1 -1], Ts);
end

% Orbit disturbance weighting matrix
Wd = ss(Wd);

% Performance weighting matrix
Wz = ss(Wz);

% Noise weighting matrix
Wn = ss(Wn);

% Port names
G.InputName = 'u';
G.OutputName = 'y';
C.InputName = 'e';
C.OutputName = 'u';
Wd.InputName = 'd';
Wd.OutputName = 'dd';
Wz.InputName = 'yd';
Wz.OutputName = 'z';
Wn.InputName = 'n';
Wn.OutputName = 'yn';

% Sum points
sum_yd = sumblk('yd = y + dd', nbpms);
sum_e = sumblk('e = r - ydn', nbpms);
sum_ydn = sumblk('ydn = yd + yn', nbpms);

P = connect(G, C, Wd, Wz, Wn, sum_yd, sum_e, sum_ydn, {'r','d','n'}, {'y','yd','z'}, {'u'});

function [sys_matrix, Ts] = sys_array2matrix(sys_array, nsys)

if ~iscell(sys_array)
    sys_array = {sys_array};
end

if length(sys_array) == 1
    corr_sys_ = cell(1, nsys);
    corr_sys_(:) = sys_array;
    sys_array = corr_sys_;
end

if length(sys_array) ~= nsys
    error('Number of elements of ''corr_sys'' must equal ''ncorr''.');
end

corr_num = cell(nsys);
corr_num(:,:) = {0};
corr_den = cell(nsys);
corr_den(:,:) = {1};
Ts_array = zeros(nsys, 1);
for i=1:nsys
    corr_tf_i = tf(sys_array{i});
    corr_tf_i = absorbDelay(corr_tf_i);
    corr_num{i,i} = corr_tf_i.num{1};
    corr_den{i,i} = corr_tf_i.den{1};
    Ts_array(i) = corr_tf_i.Ts;
end
Ts = Ts_array(1);

if ~all(Ts_array == Ts)
    error('Sample time ''Ts'' of all ''corr_sys'' systems must be identical.');
end

sys_matrix = tf(corr_num, corr_den, Ts);