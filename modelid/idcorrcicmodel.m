function [A,B,C,D,K,x0] = idcorrcicmodel(par,T,aux)

nk = aux(1);
n = aux(2);
nu = aux(3);

ncic = 2;

i=0;
for k=1:nu
    p{1+3*(k-1)} = par(i+(1:n));      i=i+length(p{end});
    p{2+3*(k-1)} = par(i+(1:n));      i=i+length(p{end});
    p{3+3*(k-1)} = par(i+1);          i=i+length(p{end});
end

Adly = [zeros(nk-1,1) eye(nk-1); 0 zeros(1,nk-1)];
Bdly = [zeros(nk-1,1); 1];
Cdly = [1 zeros(1, nk-1)];
Ddly = 0;

Acic = [0 1; 0 0];
Bcic = [0; 1];
Ccic = [0.0581814125297673 0.767679269517677];
Dcic = 1-sum(Ccic); % Enforces unitary gain

if nk ~= 0
    [Adlycic,Bdlycic,Cdlycic,Ddlycic] = concat(Adly, Bdly, Cdly, Ddly, Acic, Bcic, Ccic, Dcic);
else
    Adlycic = Acic;
    Bdlycic = Bcic;
    Cdlycic = Ccic;
    Ddlycic = Dcic;
end

for k=1:nu
    Acorr = [zeros(n-1,1) eye(n-1); p{1+3*(k-1)}'];
    Bcorr = [zeros(n-1,1); 1];
    Ccorr = p{2+3*(k-1)}';
    Dcorr = p{3+3*(k-1)};
    [A_,B_,C_,D_] = concat(Adlycic,Bdlycic,Cdlycic,Ddlycic, Acorr, Bcorr, Ccorr, Dcorr);
    A_array{k} = A_;
    B_array{k} = B_;
    C_array{k} = C_;
    D_array{k} = D_;
    K_array{k} = zeros(nk+n+ncic,1);
end

A = blkdiag(A_array{:});
B = blkdiag(B_array{:});
C = blkdiag(C_array{:});
D = blkdiag(D_array{:});
K = blkdiag(K_array{:});

x0 = zeros(nu*(n+nk+ncic),1);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, zeros(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;