clear

Ts_fofb = 10e-6;
Ts_sofb = 1e-3;
Kis = 100;
Kif = 10000;
npts = 20000;


[b,a] = butter(1, 2*100/100e3);
[b2,a2] = butter(3, 2*1/100e3);

rng(10000)

delay_fofb_bpm = 1;
delay_fofb_corr = 1;
delay_sofb_bpm = 1;
delay_sofb_corr = 0;

load sirius_matrices_2018-10-18

plane = 'y';

if strcmpi(plane, 'x')
    idx_bpm = 37 + (1:40);
    idx_fofb = 16 + (1:16);
    idx_sofb = 24 + (1:24);
    Mfofb = srius_matrices.bpm.fofb.Mx(idx_bpm, idx_fofb);
    Msofb = srius_matrices.bpm.sofb.Mx(idx_bpm, idx_sofb);
    Mdist = srius_matrices.bpm.dist.misalign.quad.Mx(idx_bpm, :);
elseif strcmpi(plane, 'y')
    idx_bpm = 196 + 37 + (1:40);
    idx_fofb = 16 + (1:16);
    idx_sofb = 32 + (1:32);
    
%     idx_bpm = 196 + 37 + (1:10);
%     idx_fofb = 16 + (1:4);
%     idx_sofb = 32 + (1:6);
    
    Mfofb = srius_matrices.bpm.fofb.My(idx_bpm, idx_fofb);
    Msofb = srius_matrices.bpm.sofb.My(idx_bpm, idx_sofb);
    Mdist = srius_matrices.bpm.dist.misalign.quad.My(idx_bpm, :);
    
    
    
    
%     Mfofb = randn(3,1);
%     Msofb = randn(3,2);
%     Mdist = randn(3,1);
else
    error('Plane not specified');
end

nbpm = size(Mfofb, 1);
ncorr_fofb = size(Mfofb, 2);
ncorr_sofb = size(Msofb, 2);
ndist = size(Mdist, 2);

Mfofb_corr = pinv(Mfofb);
Msofb_corr = pinv(Msofb);

% Analog system
sys_beam_slow = ss([], [], [], Msofb, 'inputname', 'ui_s', 'outputname', 'y_s');
sys_beam_fast = ss([], [], [], Mfofb, 'inputname', 'ui_f', 'outputname', 'y_f');
D = ss([], [], [], Mdist, 'inputname', 'd', 'outputname', 'yd');

for i=1:ncorr_sofb
    sys_corr_slow{i} = tf(1,[0.01 1], 'inputname', sprintf('u_s(%d)',i), 'outputname', sprintf('ui_s(%d)',i));
end
if ncorr_sofb == 1
    sys_corr_slow{i}.inputname = 'u_s';
    sys_corr_slow{i}.outputname = 'ui_s';
end
    
for i=1:ncorr_fofb
    sys_corr_fast{i} = tf(1,[0.0001 1], 'inputname', sprintf('u_f(%d)',i), 'outputname', sprintf('ui_f(%d)',i));
end
if ncorr_fofb == 1
    sys_corr_fast{i}.inputname = 'u_f';
    sys_corr_fast{i}.outputname = 'ui_f';
end


Ps = connect(sys_corr_slow{:}, sys_beam_slow, {'u_s'}, {'y_s'});
Pf = connect(sys_corr_fast{:}, sys_beam_fast, {'u_f'}, {'y_f'});


% FOFB controller
aux = Kif*Mfofb_corr*Ts_fofb;
fofbctrl = ss(eye(ncorr_fofb), aux, eye(ncorr_fofb), aux, Ts_fofb, 'inputname', 'ys_f', 'outputname', 'u_f_nodelay');
for i=1:nbpm
    sys_delay_fofb_bpm{i} = tf(1, 1, Ts_fofb, 'inputdelay', delay_fofb_bpm, 'inputname', sprintf('y(%d)',i), 'outputname', sprintf('ys_f(%d)',i));
end
if nbpm == 1
    sys_delay_fofb_bpm{i}.inputname = 'y';
    sys_delay_fofb_bpm{i}.outputname = 'ys_f';
end
for i=1:ncorr_fofb
    sys_delay_fofb_corr{i} = tf(1, 1, Ts_fofb, 'inputdelay', delay_fofb_corr, 'inputname', sprintf('u_f_nodelay(%d)',i), 'outputname', sprintf('u_f(%d)',i));
end
if ncorr_fofb == 1
    sys_delay_fofb_corr{i}.inputname = 'u_f_nodelay';
    sys_delay_fofb_corr{i}.outputname = 'u_f';
end
Cf = connect(fofbctrl, sys_delay_fofb_bpm{:}, sys_delay_fofb_corr{:}, 'y','u_f');

% SOFB controller
aux = Kis*Msofb_corr*Ts_sofb;
sofbctrl = ss(eye(ncorr_sofb), aux, eye(ncorr_sofb), aux, Ts_sofb, 'inputname', 'ys_s', 'outputname', 'u_s_nodelay');
for i=1:nbpm
    sys_delay_sofb_bpm{i} = tf(1, 1, Ts_sofb, 'inputdelay', delay_sofb_bpm, 'inputname', sprintf('y(%d)',i), 'outputname', sprintf('ys_s(%d)',i));
end
if nbpm == 1
    sys_delay_sofb_bpm{i}.inputname = 'y';
    sys_delay_sofb_bpm{i}.outputname = 'ys_s';
end
for i=1:ncorr_sofb
    sys_delay_sofb_corr{i} = tf(1, 1, Ts_sofb, 'inputdelay', delay_sofb_corr, 'inputname', sprintf('u_s_nodelay(%d)',i), 'outputname', sprintf('u_s(%d)',i));
end
if ncorr_sofb == 1
    sys_delay_sofb_corr{i}.inputname = 'u_s_nodelay';
    sys_delay_sofb_corr{i}.outputname = 'u_s';
end
Cs = connect(sofbctrl, sys_delay_sofb_bpm{:}, sys_delay_sofb_corr{:}, 'y','u_s');

t = (0:npts-1)'*Ts_fofb;
orbit_dist = filter(b,a,randn(npts, ndist)) + filter(b2,a2,100*randn(npts, ndist));
%fofb_bpm_noise = zeros(npts, nbpm);
%sofb_bpm_noise = zeros(npts, nbpm);

tfinal = t(end);

simsignals = [...
    t ...
    orbit_dist ...
    ];

tic
simout_simulink = sim('ofb2', 'StopTime', num2str(tfinal));
toc

r = get(simout_simulink, 'yout');
tout = r.time;
y = r.signals(1).values;
us = r.signals(2).values;
uf = r.signals(3).values;
yd = r.signals(4).values;

%plant_out = r.signals(5).values;
% y_nodist = plant_out(:, 1:nbpm);
% ui_s = plant_out(:, nbpm + (1:ncorr_sofb));
% ui_f = plant_out(:, nbpm + ncorr_sofb + (1:ncorr_fofb));
% yi_s = plant_out(:, nbpm + ncorr_sofb + ncorr_fofb + (1:nbpm));
% yi_f = plant_out(:, nbpm + ncorr_sofb + ncorr_fofb + nbpm  + (1:nbpm));

plot(tout, y, 'b'); hold all;