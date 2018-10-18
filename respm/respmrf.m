function [M, M_] = respmrf(orbit_points)

global THERING;

delta_freq = 1;

original_cavity_status = getcavity;
setcavity('on');
setradiation('on');

rf_indexes = findcells(THERING, 'Frequency');

n_orbit_points = length(orbit_points);
n_elements = length(rf_indexes);

M3 = zeros(n_orbit_points, n_elements, 4); % 3-D Matrix

for i=1:n_elements
    orbit1 = findorbit6(THERING, orbit_points);
    THERING{rf_indexes(i)}.Frequency = THERING{rf_indexes(i)}.Frequency + delta_freq/2;
    orbit2 = findorbit6(THERING, orbit_points);
    THERING{rf_indexes(i)}.Frequency = THERING{rf_indexes(i)}.Frequency - delta_freq/2;
    M3(:,i,:)  = (orbit2(1:4,:) - orbit1(1:4,:))'/delta_freq;
end

setradiation('off');
setcavity(original_cavity_status);

M = [M3(:,:,1); M3(:,:,3)];
M_ = [M3(:,:,2); M3(:,:,4)];