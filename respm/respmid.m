function [M, M_] = respmid(orbit_points, indexes, plane)

global THERING;

nss = length(indexes);

original_PassMethod = cell(nss,1);
for i=1:nss
    original_PassMethod{i} = THERING{indexes(i)}.PassMethod;
    THERING{indexes(i)}.PassMethod = 'CorrectorPass';
    THERING{indexes(i)}.KickAngle(1) = 0;
    THERING{indexes(i)}.KickAngle(2) = 0;
end

[M, M_] = respmcorr(THERING, orbit_points, indexes, plane);

for i=1:nss
    THERING{indexes(i)}.PassMethod = original_PassMethod;
end
