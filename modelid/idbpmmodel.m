function [A,B,C,D,K,x0] = idbpmmodel(par,T,aux)

nk = aux(1);
n = aux(2);
nu = aux(3);

p1 = par(1:n);
p2 = par(n+1:2*n);
p3 = par(2*n+1);

i=0;
for k=1:nu
    p{1+1*(k-1)} = par(2*n+1+k);
end

if nk ~= 0
    Adly = [zeros(nk-1,1) eye(nk-1); 0 zeros(1,nk-1)];
    Bdly = [zeros(nk-1,1); 1];
    Cdly = [1 zeros(1, nk-1)];
    Ddly = 0;
else
    Adly = [];
    Bdly = zeros(0,1);
    Cdly = zeros(1,0);
    Ddly = 1;
end

for k=1:nu
    if n ~= 0
        Abpm = [zeros(n-1,1) eye(n-1); p1'];
        Bbpm = [zeros(n-1,1); 1]*p{1+1*(k-1)};
        Cbpm = p2';
        Dbpm = p3*p{1+1*(k-1)};
    else
        Abpm = [];
        Bbpm = zeros(0,1);
        Cbpm = zeros(1,0);
        Dbpm = 1;
    end
    [A_,B_,C_,D_] = concat(Adly, Bdly, Cdly, Ddly, Abpm, Bbpm, Cbpm, Dbpm);
    A_array{k} = A_;
    B_array{k} = B_;
    C_array{k} = C_;
    D_array{k} = D_;
    K_array{k} = zeros(nk+n,1);
end

A = blkdiag(A_array{:});
B = blkdiag(B_array{:});
C = blkdiag(C_array{:});
D = blkdiag(D_array{:});
K = blkdiag(K_array{:});

x0 = zeros(nu*(n+nk),1);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, zeros(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;