function mdl = vacchambmdl(b,d,sigmac,np)
% mdl = VACCHAMBMDL(b,d,sigmac,np)
%
% b:      chamber radius [m]
% d:      chamber thickness [m]
% sigmac: chamber electrical conductivity [1/Ohm/m]
% np:     nunmber of poles of the model
%
% mdl:    vacuum chamber transfer function
%
% Based on "Eddy Current Shielding by Electrically Thick Vacuum Chambers"
% Boris Podobedov et al. (PAC'09)

mu0sigmac = 4*pi*1e-7*sigmac;
p = zeros(1,np);
p(1) = -2/(mu0sigmac*b*d);
for n=1:np-1
    p(n+1) = -(n*pi/d)^2/mu0sigmac;
end
mdl = zpk([],p,prod(abs(p)));
