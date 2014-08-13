function [A,B,C,D,K,x0] = idfofbmodel(par,T,aux)

[nbpm, ncorr, nabpm, nbbpm, nkbpm, nacorr, nbcorr, nkcorr, nacorrfix, nbcorrfix, bpmungain, corrungain, bpmhvsymmetry, M, a1n_corrfix, b0n_corrfix] = getauxvars(aux);

if bpmhvsymmetry
    nbpm2 = nbpm/2;
else
    nbpm2 = nbpm;
end

pbpmoffset = 0;
if ~bpmungain
    pcorroffset = nbpm2*(nabpm+nbbpm);
else
    pcorroffset = nbpm2*(nabpm+nbbpm-1);
end

[Abpm,Bbpm,Cbpm,Dbpm] = ssmatrices(nbpm2,nabpm,nbbpm,nkbpm,par,pbpmoffset,bpmungain);

if bpmhvsymmetry
    Abpm = blkdiag(Abpm,Abpm);
    Bbpm = blkdiag(Bbpm,Bbpm);
    Cbpm = blkdiag(Cbpm,Cbpm);
    Dbpm = blkdiag(Dbpm,Dbpm);
end

% Orbit corrector model
[Acorr,Bcorr,Ccorr,Dcorr] = ssmatrices(ncorr,nacorr,nbcorr,nkcorr,par,pcorroffset,corrungain);
if nacorrfix
    [Acorrfix,Bcorrfix,Ccorrfix,Dcorrfix] = ssmatrices(ncorr,nacorrfix,nbcorrfix,0,[a1n_corrfix b0n_corrfix],0,false);
    [Acorr,Bcorr,Ccorr,Dcorr] = concat(Acorr,Bcorr,Ccorr,Dcorr,Acorrfix,Bcorrfix,Ccorrfix,Dcorrfix);
end


% Beam response state-space matrices
Arespm = sparse(0,0);
Brespm = sparse(0,ncorr);
Crespm = sparse(nbpm,0);
Drespm = sparse(M);

[A,B,C,D] = concat(Acorr,Bcorr,Ccorr,Dcorr,Arespm,Brespm,Crespm,Drespm);
[A,B,C,D] = concat(A,B,C,D,Abpm,Bbpm,Cbpm,Dbpm);
K = zeros(nbpm*nabpm + ncorr*(nacorr+nacorrfix), nbpm);
x0 = zeros(nbpm*nabpm + ncorr*(nacorr+nacorrfix), 1);

A = full(A);
B = full(B);
C = full(C);
D = full(D);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, sparse(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;
