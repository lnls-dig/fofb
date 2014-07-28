function [A,B,C,D,K,x0] = idcorrcicmodel(par,T,aux)

nk = aux(1);
n = aux(2);
ncic = aux(3);

i = 0;
p1 = par(i+(1:n));      i=i+length(p1);
p2 = par(i+(1:n));      i=i+length(p2);
p3 = par(i+1);          i=i+length(p3);
p4 = par(i+(1:n));      i=i+length(p4);
p5 = par(i+(1:n));      i=i+length(p5);
p6 = par(i+1);          i=i+length(p6);
p7 = par(i+(1:ncic));   i=i+length(p7);
p8 = par(i+(1:ncic));   i=i+length(p8);
p9 = par(i+1);          i=i+length(p9);
p10 = par(i+(1:n));     i=i+length(p10);
p11 = par(i+(1:n));     i=i+length(p11);
p12 = par(i+(1:ncic));  i=i+length(p12);
p13 = par(i+(1:ncic));  i=i+length(p12);

Adly = [zeros(nk-1,1) eye(nk-1); 0 zeros(1,nk-1)];
Bdly = [zeros(nk-1,1); 1];
Cdly = [1 zeros(1, nk-1)];
Ddly = 0;

Acorr1 = [zeros(n-1,1) eye(n-1); p1'];
Bcorr1 = [zeros(n-1,1); 1];
Ccorr1 = p2';
Dcorr1 = p3;

Acorr2 = [zeros(n-1,1) eye(n-1); p4'];
Bcorr2 = [zeros(n-1,1); 1];
Ccorr2 = p5';
Dcorr2 = p6;

Acic = [zeros(ncic-1,1) eye(ncic-1); p7'];
Bcic = [zeros(ncic-1,1); 1];
Ccic = p8';
Dcic = p9;

[A_,B_,C_,D_] = concat(Adly, Bdly, Cdly, Ddly, Acorr1, Bcorr1, Ccorr1, Dcorr1);
[A1,B1,C1,D1] = concat(A_,B_,C_,D_, [], [], zeros(1,0), 1);

[A_,B_,C_,D_] = concat(Adly, Bdly, Cdly, Ddly, Acorr2, Bcorr2, Ccorr2, Dcorr2);
[A2,B2,C2,D2] = concat(A_,B_,C_,D_, Acic, Bcic, Ccic, Dcic);

K1 = [zeros(nk,1); p10; p12];
K2 = [zeros(nk,1); p11; p13];

A = blkdiag(A1,A2);
B = blkdiag(B1,B2);
C = blkdiag(C1,C2);
D = blkdiag(D1,D2);
K = blkdiag(K1,K2);

x0 = zeros(2*(n+nk+ncic),1);

function [A,B,C,D] = concat(A1,B1,C1,D1,A2,B2,C2,D2)

A = [A1, zeros(size(A1,2),size(A2,1)); B2*C1 A2];
B = [B1; B2*D1];
C = [D2*C1 C2];
D = D2*D1;