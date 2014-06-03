function invR = evc(R, idx, tol)

if nargin < 3 || isempty(tol)
    tol = 0;
end

% Build auxiliary matrix (C) which selects constraints
S = zeros(length(idx), size(R,1));
for i=1:length(idx)
    S(i,idx(i)) = 1;
end
C = (S*R)';

invA = pinv(R'*R, tol);
aux1 = C'*invA;
aux2 = invA*(C/(aux1*C));
invR = aux2*S-(-invA+aux2*aux1)*R';
