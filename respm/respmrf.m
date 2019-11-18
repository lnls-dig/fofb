function M = respmrf(ring, orbit_points, rf_index)

original_cavity_status = getcavity;
setcavity('on');
setradiation('on');

M = findrespm(ring, orbit_points, rf_index, 1, 'Frequency', 1,1, 'findorbit6');

setradiation('off'); % FIXME: should set radiation to previous state heve. However, there's no "getradiation" function yet to retrieve current state
setcavity(original_cavity_status);