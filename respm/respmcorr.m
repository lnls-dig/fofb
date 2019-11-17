function M = respmcorr(THERING, orbit_points, indexes, plane)

delta_exc = 1e-5;

n_orbit_points = length(orbit_points);

nelements = size(indexes,1);
M = zeros(n_orbit_points, nelements, 4);
for i=1:nelements
    if strcmpi(plane{i}, 'h')
        planenum = 1;
        planecorr = 'CH';
    else
        planenum = 2;
        planecorr = 'CV';
    end
    
    original_kick = lnls_get_kickangle(THERING, indexes, plane);
    THERING = lnls_set_kickangle(THERING, original_kick(i)-delta_exc/2, indexes(i), plane);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING = lnls_set_kickangle(THERING, original_kick(i)+delta_exc/2, indexes(i), plane);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    M(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
    THERING = lnls_set_kickangle(THERING, original_kick(i), indexes(i), plane);
end
