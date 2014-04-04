function M = fofb_orbit_respm(THERING, markers)

% Response matrix: orbit vs. corrector kicks
if ~isempty(markers.hcm) && ~isempty(markers.vcm)
    fprintf('   Response matrix: orbit vs. corrector magnet kicks\n'); tic;
    [M_bpm_hcm, M_bpm_vcm] = corrector_respm(THERING, markers.bpm, markers.hcm, markers.vcm);
    [M_source_hcm, M_source_vcm] = corrector_respm(THERING, sort(markers.source), markers.hcm, markers.vcm);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_hcm = [];
    M_bpm_vcm = [];
    M_source_hcm = [];
    M_source_vcm = [];
end

% Response matrix: orbit vs. RF frequency
if ~isempty(markers.rf)
    fprintf('   Response matrix: orbit vs. RF frequency\n'); tic;
    M_bpm_rf = rf_respm(THERING, markers.bpm, markers.rf);
    M_source_rf = rf_respm(THERING, markers.source, markers.rf);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_rf = [];
    M_source_rf = [];
end

% Response matrix: orbit vs. ID disturbance
if ~isempty(markers.id)
    fprintf('   Response matrix: orbit vs. ID disturbance\n'); tic;

    n_id = length(markers.id);
    
    original_PassMethod = cell(n_id,1);
    for i=1:n_id
        original_PassMethod{i} = THERING{markers.id(i)}.PassMethod;
        THERING{markers.id(i)}.PassMethod = 'CorrectorPass';
        THERING{markers.id(i)}.KickAngle(1) = 0;
        THERING{markers.id(i)}.KickAngle(2) = 0;
    end
    
    [M_bpm_idh, M_bpm_idv] = corrector_respm(THERING, markers.bpm, markers.id, markers.id);
    [M_source_idh, M_source_idv] = corrector_respm(THERING, markers.source, markers.id, markers.id);
    
    for i=1:n_id
        THERING{markers.id(i)}.PassMethod = original_PassMethod;
    end
            
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_idh = [];
    M_bpm_idv = [];
    M_source_idh = [];
    M_source_idv = [];
end

% Response matrix: orbit vs. sextupole displacements
if ~isempty(markers.sext)
    fprintf('   Response matrix: orbit vs. sextupoles displacements\n'); tic;
    [M_bpm_sexth, M_bpm_sextv] = displacement_respm(THERING, markers.bpm, markers.sext);
    [M_source_sexth, M_source_sextv] = displacement_respm(THERING, markers.source, markers.sext);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_sexth = [];
    M_bpm_sextv = [];
    M_source_sexth = [];
    M_source_sextv = [];
end

% Response matrix: orbit vs. quadrupole displacements
if ~isempty(markers.quad)
    fprintf('   Response matrix: orbit vs. quadrupoles displacements\n'); tic;
    [M_bpm_quadh, M_bpm_quadv] = displacement_respm(THERING, markers.bpm, markers.quad);
    [M_source_quadh, M_source_quadv] = displacement_respm(THERING, markers.source, markers.quad);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_quadh = [];
    M_bpm_quadv = [];
    M_source_quadh = [];
    M_source_quadv = [];
end

% Response matrix: orbit vs. dipole displacements
if ~isempty(markers.dipole)
    fprintf('   Response matrix: orbit vs. dipoles displacements\n'); tic;
    [M_bpm_dipoleh, M_bpm_dipolev] = displacement_respm(THERING, markers.bpm, markers.dipole);
    [M_source_dipoleh, M_source_dipolev] = displacement_respm(THERING, markers.source, markers.dipole);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    M_bpm_dipoleh = [];
    M_bpm_dipolev = [];
    M_source_dipoleh = [];
    M_source_dipolev = [];
end

% Assignment of outputs
M_bpm = struct('hcm', M_bpm_hcm, ...
               'vcm', M_bpm_vcm, ...
               'idh', M_bpm_idh, ...
               'idv', M_bpm_idv, ...
               'sexth', M_bpm_sexth, ...
               'sextv', M_bpm_sextv, ...
               'quadh', M_bpm_quadh, ...
               'quadv', M_bpm_quadv, ...
               'dipoleh', M_bpm_dipoleh, ...
               'dipolev', M_bpm_dipolev, ...
               'rf', M_bpm_rf);

M_source = struct('hcm', M_source_hcm, ...
                  'vcm', M_source_vcm, ...
                  'idh', M_source_idh, ...
                  'idv', M_source_idv, ...
                  'sexth', M_source_sexth, ...
                  'sextv', M_source_sextv, ...
                  'quadh', M_source_quadh, ...
                  'quadv', M_source_quadv, ...
                  'dipoleh', M_source_dipoleh, ...
                  'dipolev', M_source_dipolev, ...
                  'rf', M_source_rf);           
           
M = struct('bpm', M_bpm, ...
           'source', M_source);