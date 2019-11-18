function M = respmcorr(ring, orbit_indexes, corr_indexes, plane, delta_kick)

if nargin < 5 || isempty(delta_kick)
    delta_kick = 1e-6;
end

% TODO: should test the pass method of all elements, not only the element
% in the first index
if strcmp(ring{corr_indexes(1)}.PassMethod, 'CorrectorPass')
    if strcmp(plane, 'x')
        M = findrespm(ring, orbit_indexes, corr_indexes, delta_kick, 'KickAngle', 1, 1, 'findorbit4');
    elseif strcmp(plane, 'y')
        M = findrespm(ring, orbit_indexes, corr_indexes, delta_kick, 'KickAngle', 1, 2, 'findorbit4');
    end
elseif any(strcmp(ring{corr_indexes(1)}.PassMethod, {'BndMPoleSymplectic4Pass', 'BndMPoleSymplectic4RadPass', ...
        'StrMPoleSymplectic4Pass', 'StrMPoleSymplectic4RadPass', 'ThinMPolePass'}))
    if strcmp(plane, 'x')
        M = findrespm(ring, orbit_indexes, corr_indexes, -delta_kick, 'PolynomB', 1, 1, 'findorbit4');
    elseif strcmp(plane, 'y')
        M = findrespm(ring, orbit_indexes, corr_indexes, delta_kick, 'PolynomA', 1, 1, 'findorbit4');
    end
    len = getcellstruct(ring, 'Length', corr_indexes, 1,1);
end

if any(strcmp(ring{corr_indexes(1)}.PassMethod, {'CorrectorPass', 'ThinMPolePass'}))
    for i=1:length(M)
        M{i} = M{i}/delta_kick;
    end
else
    for i=1:length(M)
        M{i} = M{i}/delta_kick./repmat(len(:)', size(M{1},1), 1);
    end
end