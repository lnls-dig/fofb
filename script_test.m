if ~exist('THERING','var')
    global THERING;
    sirius;
    setoperationalmode(1);
end

displacement = 50e-9; % [m]
deltakick = 1e-4; % [rad]
deltafreq = 1000; % [Hz]

M = fofb_orbit_respm_sirius(THERING);
