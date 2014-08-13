function [nbpm, ncorr, nabpm, nbbpm, nkbpm, nacorr, nbcorr, nkcorr, nacorrfix, nbcorrfix, bpmungain, corrungain, bpmhvsymmetry, M, a1n_corrfix, b0n_corrfix] = getauxvars(auxvars)

nbpm = auxvars(1);
ncorr = auxvars(2);
nabpm = auxvars(3);
nbbpm = auxvars(4);
nkbpm = auxvars(5);
nacorr = auxvars(6);
nbcorr = auxvars(7);
nkcorr = auxvars(8);
nacorrfix = auxvars(9);
nbcorrfix = auxvars(10);
bpmungain = auxvars(11);
corrungain = auxvars(12);
bpmhvsymmetry = auxvars(13);                            i = 13;
M = reshape(auxvars(i+(1:nbpm*ncorr)), nbpm, ncorr);    i = i + nbpm*ncorr; % Static beam response parameters (beam response matrix)
a1n_corrfix = auxvars(i+(1:ncorr*nacorrfix));           i = i + ncorr*nacorrfix;
if nbcorrfix > 0
    b0n_corrfix = auxvars(i+(1:ncorr*(nbcorrfix+1)));
else
    b0n_corrfix = [];
end