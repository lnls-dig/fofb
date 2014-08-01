function [A,B,C,D,K,x0] = idfofbmodel(par,T,aux)
nbpm = aux(1);
ncorr = aux(2);
nxbpm = aux(3);
nxcorr = aux(4);
bpmungain = aux(5);
corrungain = aux(6);
M = reshape(aux(6+(1:nbpm*ncorr)), nbpm, ncorr); % Static beam response parameters (beam response matrix)

pbpmoffset = 0;
if ~bpmungain
    pcorroffset = nbpm*(2*nxbpm+1);
else
    pcorroffset = nbpm*(2*nxbpm);
end

[Abpm,Bbpm,Cbpm,Dbpm] = ssmatrices(nbpm,nxbpm,par,pbpmoffset,bpmungain);
[Acorr,Bcorr,Ccorr,Dcorr] = ssmatrices(ncorr,nxcorr,par,pcorroffset,corrungain);

% Beam response state-space matrices
Arespm = sparse(0,0);
Brespm = sparse(0,ncorr);
Crespm = sparse(nbpm,0);
Drespm = sparse(M);

[A,B,C,D] = concat(Acorr,Bcorr,Ccorr,Dcorr,Arespm,Brespm,Crespm,Drespm);
[A,B,C,D] = concat(A,B,C,D,Abpm,Bbpm,Cbpm,Dbpm);
K = zeros(nbpm*nxbpm + ncorr*(nxcorr), nbpm);
x0 = zeros(nbpm*nxbpm + ncorr*(nxcorr), 1);

A = full(A);
B = full(B);
C = full(C);
D = full(D);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, sparse(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;

function [A,B,C,D] = ssmatrices(nu,nx,par,poffset,ungain)
nux = nu*nx;
if nx ~= 0
    aux1 = reshape(repmat(0:nx:nux-1,nx,1), 1, []);
    A = sparse([setdiff(1:nux, 1:nx:nux) 1:nux], ...
        [setdiff(1:nux, nx:nx:nux) repmat(nx,1,nux)+aux1], ...
        [ones(1,nu*(nx-1)) par(poffset+(1:nux))']);
    B = sparse([1:nx:nux nux], [1:nu nu], [ones(1,nu) 0]);
    C = sparse(reshape(repmat(1:nu,nx,1),1,[]), 1:nux, par(poffset+(nux+1:2*nux))');
else
    A = sparse(0,0);
    B = sparse(0,nu);
    C = sparse(nu,0);
end
if ~ungain
    D = sparse(1:nu, 1:nu, par(poffset+(2*nux+1:2*nux+nu))');
else
    D = speye(nu)-C*((speye(nux)-A)\B);
end
