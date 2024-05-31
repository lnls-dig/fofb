function T = ofbmdl(M, Mc, K, corr_sys, F, Wd, Wz)

%TODO: add shaping filters - F(z)

nbpms = size(M,1);
ncorr = size(M,2);

if nargin < 6 || isempty(Wd)
    Wd = eye(nbpms);
end

if nargin < 7 || isempty(Wz)
    Wz = eye(nbpms);
end

% Actuator transfer function (includes power supply, magent, vacuum
% chamber, BPM, beam and network delay dynamics)
ncorr_sys = length(corr_sys);
corr_num = cell(ncorr);
corr_num(:,:) = {0};
corr_den = cell(ncorr);
corr_den(:,:) = {1};
for i=1:ncorr_sys
    corr_tf_i = tf(corr_sys{i});
    corr_tf_i = absorbDelay(corr_tf_i);
    corr_num{i,i} = corr_tf_i.num{1};
    corr_den{i,i} = corr_tf_i.den{1};
    %corr_delay(i,i) = totaldelay(corr_tf_i);
end
Ts = corr_tf_i.Ts;
%corr_tf = tf(corr_num, corr_den, Ts, 'ioDelay', corr_delay);
corr_tf = tf(corr_num, corr_den, Ts);

% Plant transfer function
G = ss(M)*ss(corr_tf);

% Controller transfer function
C = diag(K)*Mc*tf([Ts 0], [1 -1], Ts);

% Orbit disturbance transfer function
Wd = ss(Wd);

% Performance transfer function
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

T = connect(G, C, Wd, Wz, sum_yd, sum_e, {'r','d'}, {'yd','z'}, {'u'});