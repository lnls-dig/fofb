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
    
    if isfield(THERING{indexes(i,:)}, 'KickAngle')
        original_kick = THERING{indexes(i,:)}.KickAngle(planenum);
        THERING{indexes(i,:)}.KickAngle(planenum) = original_kick - (delta_exc/2);
        orbit1 = findorbit4(THERING, 0, orbit_points);
        THERING{indexes(i,:)}.KickAngle(planenum) = original_kick + (delta_exc/2);
        orbit2 = findorbit4(THERING, 0, orbit_points);
        THERING{indexes(i,:)}.KickAngle(planenum) = original_kick;
    elseif isfield(THERING{indexes(i,:)}, 'CH') || isfield(THERING{indexes(i,:)}, 'CV')
        original_kick = THERING{indexes(i,:)}.(planecorr);
        THERING{indexes(i,:)}.(planecorr) = original_kick - (delta_exc/2);
        orbit1 = findorbit4(THERING, 0, orbit_points);
        THERING{indexes(i,:)}.(planecorr) = original_kick + (delta_exc/2);
        orbit2 = findorbit4(THERING, 0, orbit_points);
        THERING{indexes(i,:)}.(planecorr) = original_kick;
    end

    M(:,i,:)  = (orbit2 - orbit1)'/delta_exc;
end
