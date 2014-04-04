function [Mx, My] = corrector_respm(THERING, orbit_points, indexes_x, indexes_y)

delta_exc = 1e-5;

n_orbit_points = length(orbit_points);

n_elements_x = size(indexes_x,1);
Mx = zeros(n_orbit_points, n_elements_x, 4);
for i=1:n_elements_x
    original_kick = THERING{indexes_x(i,:)}.KickAngle(1);
    THERING{indexes_x(i,:)}.KickAngle(1) = original_kick - (delta_exc/2);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{indexes_x(i,:)}.KickAngle(1) = original_kick + (delta_exc/2);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    THERING{indexes_x(i,:)}.KickAngle(1) = original_kick;

    Mx(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
end

n_elements_y = size(indexes_y,1);
My = zeros(n_orbit_points, n_elements_y, 4);
for i=1:n_elements_y
    original_kick = THERING{indexes_y(i,:)}.KickAngle(2);
    THERING{indexes_y(i,:)}.KickAngle(2) = original_kick - (delta_exc/2);
    orbit1 = findorbit4(THERING, 0, orbit_points);
    THERING{indexes_y(i,:)}.KickAngle(2) = original_kick + (delta_exc/2);
    orbit2 = findorbit4(THERING, 0, orbit_points);
    THERING{indexes_y(i,:)}.KickAngle(2) = original_kick;
    
    My(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
end