function [sys_bpm, sys_corr] = idparam2tf(param, nk, nbpm, ncorr, nxbpm, nxcorr, bpmungain, corrungain)

pbpmoffset = 0;
pcorroffset = 2*nxbpm+1;

p1bpm = param(pbpmoffset + (1:nbpm*nxbpm));
p2bpm = param(pbpmoffset + (nbpm*nxbpm+1:nbpm*2*nxbpm));
if ~bpmungain
    p3bpm = param(pbpmoffset + (nbpm*2*nxbpm+1:nbpm*(2*nxbpm+1)));
else
    p3bpm = [];
end
sys_bpm = idparam2tf_(p1bpm, p2bpm, p3bpm, 0, nbpm, nxbpm, bpmungain);

p1corr = param(pcorroffset + (1:ncorr*nxcorr));
p2corr = param(pcorroffset + (ncorr*nxcorr+1:ncorr*2*nxcorr));
if ~corrungain
    p3corr = param(pcorroffset + (ncorr*2*nxcorr+1:ncorr*(2*nxcorr+1)));
else
    p3corr = [];
end
sys_corr = idparam2tf_(p1corr, p2corr, p3corr, nk, ncorr, nxcorr, corrungain);

function sys = idparam2tf_(p1, p2, p3, nk, nu, nx, ungain)

Ts = 320e-6;

for i=1:nu
    A = [[zeros(1,nx-1); eye(nx-1)] p1((i-1)*nx + (1:nx))];
    if nx > 0
        B = [1; zeros(nx-1,1)];
    else
        B = zeros(0,1);
    end
    C = p2((i-1)*nx + (1:nx))';
    if ~ungain
        D = p3((i-1) + 1);
    else
        D = 1-C*((eye(nx)-A)\B);
    end
    sys{i} = ss(A,B,C,D,Ts);
    if nk > 0
        sys{i}.InputDelay = nk;
    end
end