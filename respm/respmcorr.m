function M = respmcorr(THERING, orbit_points, indexes, plane, delta_exc)

if nargin < 5 || isempty(delta_exc)
    delta_exc = 1e-5;
end

n_orbit_points = length(orbit_points);

nelements = size(indexes,1);
M = zeros(n_orbit_points, nelements, 4);

original_kick = lnls_get_kickangle(THERING, indexes, plane);
for i=1:nelements
    kicks = original_kick;
    THERING = lnls_set_kickangle(THERING, original_kick(i)-delta_exc/2, indexes(i), plane);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    kicks(i) = kicks(i)-delta_exc/2;
    THERING = lnls_set_kickangle(THERING, original_kick(i)+delta_exc/2, indexes(i), plane);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    M(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
    THERING = lnls_set_kickangle(THERING, original_kick(i), indexes(i), plane);
end