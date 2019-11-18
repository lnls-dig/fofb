function M = respmid(ring, orbit_indexes, id_indexes, plane, delta_kick)

n = length(id_indexes);

if nargin < 5 || isempty(delta_kick)
    delta_kick = 1e-6;
end

original_PassMethod = cell(n,1);
for i=1:n
    original_PassMethod{i} = ring{id_indexes(i)}.PassMethod;
    ring{id_indexes(i)}.PassMethod = 'CorrectorPass';
    ring{id_indexes(i)}.KickAngle = [0 0];
end

M = respmcorr(ring, orbit_indexes, id_indexes, plane, delta_kick);

for i=1:n
    ring{id_indexes(i)}.PassMethod = original_PassMethod;
end
