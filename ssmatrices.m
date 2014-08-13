function [A,B,C,D] = ssmatrices(nu,na,nb,nk,par,poffset,ungain)
nx = max(na,nb+nk-1);

nux = nu*nx;
nub = nu*nb;
if ungain
    nfix=1;
else
    nfix=0;
end

a1n = reshape(par(poffset+(1:nu*na)), na, nu);
b0n = reshape(par(poffset+nu*na+(1:nub-nu*nfix)), nb-nfix, nu);
if nfix
    b0n = [(sum(a1n,1)+1)-sum(b0n,1); b0n];
end

a1n = [a1n; zeros(nx-na,nu)];
b0n = [zeros(nk,nu); b0n; zeros(nx+1-nk-nb,nu)];

b0 = b0n(1,:);
b1n = b0n(2:nx+1,:);

ia1n = a1n(end:-1:1,:);
ib1n = b1n(end:-1:1,:);
b0 = b0(:)';
ia1n = ia1n(:)';
ib1n = ib1n(:)';

if nx ~= 0
    aux1 = reshape(repmat(0:nx:nux-1,nx,1), 1, []);
    A = sparse([setdiff(1:nux, nx:nx:nux) repmat(nx,1,nux)+aux1], ...
        [setdiff(1:nux, 1:nx:nux) 1:nux], ...
        [ones(1,nu*(nx-1)) -ia1n],nux,nux);
    B = sparse(nx:nx:nux, 1:nu, ones(1,nu), nux, nu);
    C = sparse(reshape(repmat(1:nu,nx,1),1,[]), 1:nux, ib1n-ia1n.*reshape(repmat(b0,nx,1),1,[]), nu, nux);
    D = sparse(1:nu, 1:nu, b0);
else
    A = sparse(0,0);
    B = sparse(0,nu);
    C = sparse(nu,0);
    D = sparse(1:nu, 1:nu, b0);
end