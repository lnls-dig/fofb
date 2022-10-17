function [T,M,Mc] = fofb_sofb_tf(param, Wd, Wz)

% TODO: must check that the number of BPMs is the same in all response matrices
nbpms = size(param(1).M,1);

if nargin < 2 || isempty(Wd)
    Wd = eye(nbpms);
end

if nargin < 3 || isempty(Wz)
    Wz = eye(nbpms);
end

for i=1:length(param)
    M{i} = param(i).M;
    Mc{i} = param(i).Mc;
    
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
    C{i}.OutputName = ['u' param_char];
    H{i}.InputName = 'e';
    H{i}.OutputName = ['ee' param_char];
    D{i}.InputName = ['d' param_char];
    D{i}.OutputName = ['dd' param_char];
    
    % Sum points
    sum_yd{i} = sumblk(sprintf('yd%s = y%s + dd%s', param_char, param_char, param_char), nbpms);
end

% Orbit disturbance transfer function
Wd = ss(Wd);
Wd.InputName = 'd';
Wd.OutputName = 'dd';

% Performance transfer function
Wz = ss(Wz);
Wz.InputName = 'ydd';
Wz.OutputName = 'z';

sum_e = sumblk('e = r - ydd', nbpms);
sum_ydsf = sumblk('yd = yds + ydf', nbpms);
sum_ydd = sumblk('ydd = yd + dd', nbpms);

T = connect(G{:}, C{:}, D{:}, H{:}, Wd, Wz, sum_yd{:}, sum_e, sum_ydsf, sum_ydd, {'r','d','ds','df'}, {'ydd','z'}, {'us','uf'});
%T = minreal(T);

end
