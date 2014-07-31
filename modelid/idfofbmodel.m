function [A,B,C,D,K,x0] = idfofbmodel(par,T,aux)
nk = aux(1);
nbpm = aux(2);
ncorr = aux(3);
nxbpm = aux(4);
nxcorr = aux(5);
M = reshape(aux(5+(1:nbpm*ncorr)), nbpm, ncorr); % Static beam response parameters (beam response matrix)

pbpmoffset = 1;
pcorroffset = (2*nxbpm+1)+1;

if nk ~= 0
    aux1 = reshape(repmat(0:nk:nk*ncorr-1, nk-1,1), 1, []);
   
    Adly = sparse([(repmat(1:(nk-1), 1, ncorr)+aux1) nk*ncorr], [(repmat(2:nk, 1, ncorr)+aux1) nk*ncorr], [ones(1,ncorr*(nk-1)) 0]);
	Bdly = sparse(nk*(1:ncorr), 1:ncorr, 1);
    Cdly = sparse([1:ncorr ncorr], [1:nk:ncorr*nk ncorr*nk], [ones(1,ncorr) 0]);
    Ddly = sparse(ncorr, ncorr);
else
    Adly = sparse(0,0);
    Bdly = sparse(0,ncorr);
    Cdly = sparse(ncorr,0);
    Ddly = speye(ncorr);
end

[Abpm,Bbpm,Cbpm,Dbpm] = ssmatrices(nbpm,nxbpm,par,pbpmoffset);
[Acorr,Bcorr,Ccorr,Dcorr] = ssmatrices(ncorr,nxcorr,par,pcorroffset);

% Beam response state-space matrices
Arespm = sparse(0,0);
Brespm = sparse(0,ncorr);
Crespm = sparse(nbpm,0);
Drespm = sparse(M);

[A,B,C,D] = concat(Acorr,Bcorr,Ccorr,Dcorr,Adly,Bdly,Cdly,Ddly);
[A,B,C,D] = concat(A,B,C,D,Arespm,Brespm,Crespm,Drespm);
[A,B,C,D] = concat(A,B,C,D,Abpm,Bbpm,Cbpm,Dbpm);
K = zeros(nbpm*nxbpm + ncorr*(nk+nxcorr), nbpm);
x0 = zeros(nbpm*nxbpm + ncorr*(nxcorr+nk), 1);

A = full(A);
B = full(B);
C = full(C);
D = full(D);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, sparse(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;

function [A,B,C,D] = ssmatrices(nu,nx,par,poffset)
if nx ~= 0
    nux = nu*nx;
    aux1 = reshape(repmat(0:nx:nux-1,nx,1), 1, []);
    A = sparse([setdiff(1:nux, 1:nx:nux) 1:nux], ...
        [setdiff(1:nux, nx:nx:nux) repmat(nx,1,nux)+aux1], ...
        [ones(1,nu*(nx-1)) par(poffset+(0:nux-1))']);
	B = sparse([1:nx:nux nux], [1:nu nu], [ones(1,nu) 0]);
    C = sparse(reshape(repmat(1:nu,nx,1),1,[]), 1:nux, par(poffset+(nux:2*nux-1))');
    D = sparse(1:nu, 1:nu, par(poffset+(2*nux:2*nux+nu-1))');
else
    A = sparse(0,0);
    B = sparse(0,nu);
    C = sparse(nu,0);
    D = speye(nu);
end