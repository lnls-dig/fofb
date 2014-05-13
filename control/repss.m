function [A,B,C,D] = repss(A,B,C,D,n)
%REPSS Replicate state-space model N times.
%
%   [A,B,C,D] = repss(A,B,C,D,n)

if ~iscell(A) && ~iscell(B) && ~iscell(C) && ~iscell(D)
    if nargin >= 5 && isscalar(n) && n > 0
        [Acell{1:n}] = deal(A);
        [Bcell{1:n}] = deal(B);
        [Ccell{1:n}] = deal(C);
        [Dcell{1:n}] = deal(D);
    else
        error('Must specify how many times matrices should be replicated.')
    end
elseif iscell(A) && iscell(B) && iscell(C) && iscell(D)
    if nargin < 5
        Acell = A;
        Bcell = B;
        Ccell = C;
        Dcell = D;
    else
        error('Cannot replicate matrices.')
    end
else
    error('A,B,C,D inputs must be all matrices or all cell of matrices.');
end

A = blkdiag(Acell{:});
B = blkdiag(Bcell{:});
C = blkdiag(Ccell{:});
D = blkdiag(Dcell{:});