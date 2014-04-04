function [Mpseudoinv, sv, newsv] = pseudoinv(M, discardfcn)
%PSEUDOINV Calculate the pseudo inverse of the matrix M using the Singular
%   Value Decomposition (SVD) method. Optionally a function DISCARDFCN can
%   be applied to the singular values in order to shape or discard the
%   singular values.
%
%   [Mpseudoinv, sv, newsv] = pseudoinv(M, discardfcn)
%   
%   Inputs
%       M: matrix
%       DISCARDFCN (optional): function handle for shaping or discarding
%           singular values. To discard a singular values a 0 must be
%           written.
%
%   Outputs
%       MPSEUDOINV: pseudo inverse matrix
%       SV: original singular values
%       NEWSV: singular values taken into account for the pesudo inverse
%           calculation.

[U,S,V] = svd(M);

invS = zeros(size(S))';
nsv = min(size(M));
sv = diag(S);
newsv = sv;

if nargin >= 2
    newsv = discardfcn(sv);
end

invS(1:nsv,1:nsv) = diag(1./newsv);
invS(invS == Inf) = 0;

Mpseudoinv = V*invS*U';