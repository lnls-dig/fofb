function [A,B,C,D,K,x0] = idfofbmodel(par,T,aux)
nbpm = aux(1);
ncorr = aux(2);
nabpm = aux(3);
nbbpm = aux(4);
nacorr = aux(5);
nbcorr = aux(6);
bpmungain = aux(7);
corrungain = aux(8);
M = reshape(aux(8+(1:nbpm*ncorr)), nbpm, ncorr); % Static beam response parameters (beam response matrix)

pbpmoffset = 0;
if ~bpmungain
    pcorroffset = nbpm*(nabpm+nbbpm+1);
else
    pcorroffset = nbpm*(nabpm+nbbpm);
end

bpmhvsymmetry = false; % TODO: include bpmhvsymmetry as paramter
if bpmhvsymmetry
    [Abpm,Bbpm,Cbpm,Dbpm] = ssmatrices(nbpm/2,nabpm,nbbpm,par,pbpmoffset,bpmungain);
    Abpm = blkdiag(Abpm,Abpm);
    Bbpm = blkdiag(Bbpm,Bbpm);
    Cbpm = blkdiag(Cbpm,Cbpm);
    Dbpm = blkdiag(Dbpm,Dbpm);
else
    [Abpm,Bbpm,Cbpm,Dbpm] = ssmatrices(nbpm,nabpm,nbbpm,par,pbpmoffset,bpmungain);
end
[Acorr,Bcorr,Ccorr,Dcorr] = ssmatrices(ncorr,nacorr,nbcorr,par,pcorroffset,corrungain);



% Beam response state-space matrices
Arespm = sparse(0,0);
Brespm = sparse(0,ncorr);
Crespm = sparse(nbpm,0);
Drespm = sparse(M);

[A,B,C,D] = concat(Acorr,Bcorr,Ccorr,Dcorr,Arespm,Brespm,Crespm,Drespm);
[A,B,C,D] = concat(A,B,C,D,Abpm,Bbpm,Cbpm,Dbpm);
K = zeros(nbpm*nabpm + ncorr*(nacorr), nbpm);
x0 = zeros(nbpm*nabpm + ncorr*(nacorr), 1);

A = full(A);
B = full(B);
C = full(C);
D = full(D);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, sparse(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;

function [A,B,C,D] = ssmatrices(nu,na,nb,par,poffset,ungain)
nua = nu*na;
nub1 = nu*(nb+1);
if ungain
    nfix=1;
else
    nfix=0;
end

a1n = reshape(par(poffset+(1:nua)), na, nu);
b0n = reshape(par(poffset+(nua+1:nua+nub1-nu*nfix)), nb+1-nfix, nu);
if nfix
    b0n = [(sum(a1n,1)+1)-sum(b0n,1); b0n];
end
b0n = [zeros(na-nb,nu); b0n];
b0 = b0n(1, :);
b1n = b0n(2:na+1,:);

ia1n = a1n(end:-1:1,:);
ib1n = b1n(end:-1:1,:);
b0 = b0(:)';
ia1n = ia1n(:)';
ib1n = ib1n(:)';

if na ~= 0
    aux1 = reshape(repmat(0:na:nua-1,na,1), 1, []);
    A = sparse([setdiff(1:nua, na:na:nua) repmat(na,1,nua)+aux1], ...
        [setdiff(1:nua, 1:na:nua) 1:nua], ...
        [ones(1,nu*(na-1)) -ia1n],nua,nua);
    B = sparse(na:na:nua, 1:nu, ones(1,nu), nua, nu);
    C = sparse(reshape(repmat(1:nu,na,1),1,[]), 1:nua, ib1n-ia1n.*reshape(repmat(b0,na,1),1,[]), nu, nua);
    D = sparse(1:nu, 1:nu, b0);
else
    A = sparse(0,0);
    B = sparse(0,nu);
    C = sparse(nu,0);
    D = sparse(1:nu, 1:nu, b0);
end