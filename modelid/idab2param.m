function param = idab2param(a,b,order,ungain)

if nargin < 4
    ungain = false;
end

na = order(1);
nb = order(2);
nk = order(3);

b = b(:,nk+(1:nb));

a0 = a(:,1);
aparam = a(:,2:na+1)./repmat(a0,1,na);
bparam = b./repmat(a0,1,nb);

if ungain
   bparam = bparam(:,2:end);
end

aparam = aparam';
bparam = bparam';

param = [aparam(:)' bparam(:)'];