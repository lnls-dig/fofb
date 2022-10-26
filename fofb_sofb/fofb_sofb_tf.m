function [T,M,Mc] = fofb_sofb_tf(param, Wd, Wz, kx)

% TODO: must check that the number of BPMs is the same in all response matrices
nbpms = size(param(1).M,1);

if nargin < 2 || isempty(Wd)
    Wd = eye(nbpms);
end

if nargin < 3 || isempty(Wz)
    Wz = eye(nbpms);
end

if nargin < 4 || isempty(kx)
    kx = 0;
end

for i=1:length(param)
    M{i} = param(i).M;
    Mc{i} = param(i).Mc;
    ncorr = size(param(i).M,2);
    
    % Transfer functions
    G{i} = M{i}*tf(1,1,'iodelay', param(i).dly);
    C{i} = tf(param(i).Ki,[1 0])*Mc{i};
    H{i} = param(i).H*eye(nbpms);
    D{i} = ss([],[],[],M{i}/norm(M{i}));

    % Limit low frequency gain of integrator to avoid numerical issues
    C{i} = C{i}*tf([1 0],[1 2*pi*1e-12]);

    param_char = param(i).char;
    
    % Port names
    G{i}.InputName = ['u' param_char];
    G{i}.OutputName = ['y' param_char];
    C{i}.InputName = ['ee' param_char];
    C{i}.OutputName = ['ui' param_char];
    H{i}.InputName = 'e';
    H{i}.OutputName = ['ee' param_char];
    D{i}.InputName = ['d' param_char];
    D{i}.OutputName = ['dd' param_char];
    
    % Sum points
    sum_yd{i} = sumblk(sprintf('yd%s = y%s + dd%s', param_char, param_char, param_char), nbpms);
    sum_ue{i} = sumblk(sprintf('u%s = ui%s + ue%s', param_char, param_char, param_char), ncorr);
end

% Orbit disturbance transfer function
Wd = ss(Wd);
Wd.InputName = 'd';
Wd.OutputName = 'dd';

% Performance transfer function
Wz = ss(Wz);
Wz.InputName = 'ydd';
Wz.OutputName = 'z';

% FOFB-to-SOFB actuation download transfer function
X = tf(kx, [1 0])*Mc{1}*M{2};
X.InputName = 'uf';
X.OutputName = 'ues';

sum_e = sumblk('e = r - ydd', nbpms);
sum_ydsf = sumblk('yd = yds + ydf', nbpms);
sum_ydd = sumblk('ydd = yd + dd', nbpms);

T = connect(G{:}, C{:}, D{:}, H{:}, Wd, Wz, X, sum_yd{:}, sum_ue{:}, sum_e, sum_ydsf, sum_ydd, {'r','d','ds','df','uef'}, {'ydd','z'}, {'us','uf'});
%T = minreal(T);

end
