function [a,b] = idparam2ab(param, nu, order, ungain)

if nargin < 4
    ungain = false;
end

na = order(1);
nb = order(2);
nk = order(3);

a = zeros(nu, max(na+1,nb+nk));
b = zeros(nu, max(na+1,nb+nk));

a(:,1:na+1) = [ones(nu,1) reshape(param(1:nu*na), na, [])'];

if ~ungain
    b(:,nk+1:nk+nb) = reshape(param(nu*na+(1:nu*nb)), nb, [])';
else
    b_ = reshape(param(nu*na+(1:nu*(nb-1))), nb-1, [])';
    b0 = sum(a,2) - sum(b_,2);
    b(:,nk+1:nk+nb) = [b0 b_];
end



