function [a,b,c,d] = repss(a,b,c,d,n)
%REPSS Replicate state-space model N times.
%
%   [a,b,c,d] = repss(a,b,c,d,n)

[acell{1:n}] = deal(a);
a = blkdiag(acell{:});

[bcell{1:n}] = deal(b);
b = blkdiag(bcell{:});

[ccell{1:n}] = deal(c);
c = blkdiag(ccell{:});

[dcell{1:n}] = deal(d);
d = blkdiag(dcell{:});