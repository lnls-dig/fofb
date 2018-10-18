function [M, M_] = respmcorr(orbit_points, indexes, plane, delta_exc)

global THERING;

if nargin < 5 || isempty(delta_exc)
    delta_exc = 1e-5;
end

n_orbit_points = length(orbit_points);

n_elements = size(indexes,1);
M3 = zeros(n_orbit_points, n_elements, 4); % 3-D Matrix

original_kick = lnls_get_kickangle(THERING, indexes, plane);
for i=1:n_elements
    THERING = lnls_set_kickangle(THERING, original_kick(i)-delta_exc/2, indexes(i), plane);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING = lnls_set_kickangle(THERING, original_kick(i)+delta_exc/2, indexes(i), plane);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    M3(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
    THERING = lnls_set_kickangle(THERING, original_kick(i), indexes(i), plane);
end

M = [M3(:,:,1); M3(:,:,3)];
M_ = [M3(:,:,2); M3(:,:,4)];