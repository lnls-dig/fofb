function [P, G, C] = ofbmdl(M, Mc, K, A, F, Wd, Wz)
% OFBMDL Orbit Feedback Model.
% 
% Construct orbit feedback state-space model.
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
%   F:    Shaping filters per orbit corrector (dynamic system object (tf,
%         zpk, ss, etc.) or cell array of dynamic system objects with NCORR
%         elements). If only one filter is provided the shaping filters
%         of all actuators are assumed to be identical.
%   Wd:   Disturbance weighting matrix, used for loop optimzation purposes.
%   Wz:   Performance weighting matrix, used for loop optimzation purposes.
%
% OUTPUTS:
%   P:    Generalized plant of the closed-loop orbit feedback system.
%   G:    Combined responses of actuators and beam (only static beam
%         response).
%   C:    Combined responses of the controller and shaping filters.

nbpms = size(M,1);
ncorr = size(M,2);

if nargin < 5 || isempty(F)
    F = tf(1, 1, -1);
end

if nargin < 6 || isempty(Wd)
    Wd = eye(nbpms);
end

if nargin < 7 || isempty(Wz)
    Wz = eye(nbpms);
end

% Actuator transfer function (includes power supply, magnet, vacuum
% chamber, BPM, beam and network delay dynamics)
[A, Ts] = sys_array2matrix(A, ncorr);

% Shaping filters
[F, Ts_F] = sys_array2matrix(F, ncorr);

if Ts_F ~= -1 && Ts ~= Ts_F
    error('Sample time ''Ts'' of actuator transfer functions and shaping filters must be equal.');
end

% Plant transfer function
G = ss(M)*ss(A);

% Controller transfer function
C = ss(F)*diag(K)*Mc*tf([Ts 0], [1 -1], Ts);

% Orbit disturbance weighting matrix
Wd = ss(Wd);

% Performance weighting matrix
Wz = ss(Wz);

% Port names
G.InputName = 'u';
G.OutputName = 'y';
C.InputName = 'e';
C.OutputName = 'u';
Wd.InputName = 'd';
Wd.OutputName = 'dd';
Wz.InputName = 'yd';
Wz.OutputName = 'z';

% Sum points
sum_yd = sumblk('yd = y + dd', nbpms);
sum_e = sumblk('e = r - yd', nbpms);

P = connect(G, C, Wd, Wz, sum_yd, sum_e, {'r','d'}, {'yd','z'}, {'u'});

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