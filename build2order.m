function [num, den] = build2order(bw, damping)

omegan = 2*pi*bw/sqrt(1 - 2*damping^2 + sqrt(2-4*damping^2+4*damping^4));
num = omegan^2;
den = [1 2*damping*omegan omegan^2];