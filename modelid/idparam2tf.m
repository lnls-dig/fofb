function [sys_bpm, sys_corr] = idparam2tf(param, modelinfo)

nbpm = modelinfo.nbpm;
ncorr = modelinfo.ncorr;
nabpm = modelinfo.nabpm;
nbbpm = modelinfo.nbbpm;
nacorr = modelinfo.nacorr;
nbcorr = modelinfo.nbcorr;
bpmungain = modelinfo.bpmungain;
corrungain = modelinfo.corrungain;
bpmhvsymmetry = modelinfo.bpmhvsymmetry;
nk = modelinfo.nk;

if bpmungain
    bpmnfix=1;
else
    bpmnfix=0;
end
if corrungain
    corrnfix=1;
else
    corrnfix=0;
end

pbpmoffset = 0;
pcorroffset = nabpm+nbbpm+1-bpmnfix;

if bpmhvsymmetry
    nbpm = nbpm/2;
end

bpm_a = [ones(1,nbpm); reshape(param(pbpmoffset + (1:nbpm*nabpm)), nabpm, [])];
bpm_b = reshape(param(pbpmoffset + (nbpm*nabpm+1:nbpm*nabpm+nbpm*(nbbpm+1)-nbpm*bpmnfix)), nbbpm+1-bpmnfix, []);
if isempty(bpm_b)
    bpm_b = zeros(0,nbpm);
end
if bpmungain
    bpm_blast = sum(bpm_a,1)-sum(bpm_b,1);
else
    bpm_blast = [];
end
bpm_b = [bpm_blast; bpm_b];
sys_bpm = idparam2tf_(bpm_a, bpm_b, 0, nbpm);

corr_a = [ones(1,ncorr); reshape(param(pcorroffset + (1:ncorr*nacorr)), nacorr, [])];
corr_b = reshape(param(pcorroffset + (ncorr*nacorr+1:ncorr*nacorr+ncorr*(nbcorr+1)-ncorr*corrnfix)), nbcorr+1-corrnfix, []);
if isempty(corr_b)
    corr_b = zeros(0,ncorr);
end
if corrungain
    corr_blast = sum(corr_a,1)-sum(corr_b,1);
else
    corr_blast = [];
end
corr_b = [corr_blast; corr_b];
sys_corr = idparam2tf_(corr_a, corr_b, nk, ncorr);

function sys = idparam2tf_(a, b, nk, nu)

Ts = 320e-6;

for i=1:nu
    sys{i} = tf(b(:,i)',a(:,i)',Ts);
    if nk > 0
        sys{i}.InputDelay = nk;
    end
end