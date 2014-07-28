function [A,B,C,D,K,x0] = iduvxmodel(par,T)

Acorr_cic
Bcorr_cic
Ccorr_cic
Dcorr_cic

Amagf
Bmagf
Cmagf
Dmagf

Abeam
Bbeam
Cbeam
Dbeam

Abpm
Bbpm
Cbpm
Dbpm



A = [0 1; 0 par(1)]; 
B = [0;par(2)];
C = eye(2);
D = zeros(2,1);
K = zeros(2,2);
x0 =[par(3);0];