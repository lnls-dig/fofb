function M = respmrf(ring, orbit_points, rf_index)

dp = 1e-9;
alphac = mcf(ring);

freq = ring{rf_index}.Frequency;

orbit = findorbit4(ring,dp,1:length(ring)+1);
M_array = -orbit(:,orbit_points)'/dp/alphac/freq;
M = {M_array(:,1) M_array(:,2) M_array(:,3) M_array(:,4)};