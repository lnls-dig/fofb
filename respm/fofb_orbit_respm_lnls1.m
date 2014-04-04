function M = fofb_orbit_respm_lnls1(THERING)

% Markers
bpm_markers = findcells(THERING, 'FamName', 'BPM')';

rad4g_markers = findcells(THERING, 'FamName', 'RAD4G')';
rad15g_markers = findcells(THERING, 'FamName', 'RAD15G')';
id_markers = findcells(THERING, 'FamName', 'LCENTER')';
id_markers = id_markers([1 4:end]); % Remove injection and cavity sections
source_markers =  [rad4g_markers; rad15g_markers; id_markers];

hcm_markers = findcells(THERING, 'FamName', 'HCM')';
vcm_markers = findcells(THERING, 'FamName', 'VCM')';
quad_markers = [];
dipole_markers = [];
rf_markers = findcells(THERING, 'FamName', 'RF')';

markers = struct('bpm', bpm_markers, ...
                 'source', source_markers, ...
                 'hcm', hcm_markers, ...
                 'vcm', vcm_markers, ...
                 'id', id_markers, ...
                 'quad', quad_markers, ...
                 'dipole', dipole_markers, ...
                 'rf', rf_markers);

M = fofb_orbit_respm(THERING, markers);

