function [Mx, My] = displacement_respm(THERING, orbit, indexes)

delta_mis = 1e-6;

n_orbit_points = length(orbit);
n_elements = size(indexes,1);

Mx = zeros(n_orbit_points, n_elements, 4);
My = zeros(n_orbit_points, n_elements, 4);

TR0 = THERING;
for i=1:n_elements
    idx = indexes(i,:);
    
    THERING = TR0;
    THERING = lnls_set_misalignmentX(ones(size(idx))*delta_mis/2, idx, THERING);
    orbit1 = findorbit4(THERING, 0, orbit);
    THERING = lnls_set_misalignmentX(-ones(size(idx))*delta_mis/2, idx, THERING);
    orbit2 = findorbit4(THERING, 0, orbit);
    Mx(:,i,:) = (orbit2-orbit1)'/delta_mis;
    
    THERING = TR0;
    THERING = lnls_set_misalignmentY(ones(size(idx))*delta_mis/2, idx, THERING);
    orbit1 = findorbit4(THERING, 0, orbit);
    THERING = lnls_set_misalignmentY(-ones(size(idx))*delta_mis/2, idx, THERING);
    orbit2 = findorbit4(THERING, 0, orbit);
    My(:,i,:) = (orbit2-orbit1)'/delta_mis;
end