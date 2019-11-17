function [Mcorr, Mss, Mdisp, Mrf] = respmorbit(THERING, markers, plane)

if strcmpi(plane, 'x')
    index_pos = 1;
    index_ang = 2;
elseif strcmpi(plane, 'y')
    index_pos = 3;
    index_ang = 4;
else
    error('Value assigned to plane is not valid');
end

% Response matrix: orbit vs. corrector kicks
if isfield(markers, 'corr') && ~isempty(markers.corr)
    fprintf('   Response matrix: orbit vs. corrector magnet kicks\n'); tic;
    Mcorr_ = respmcorr(THERING, markers.orbit, markers.corr, plane);
    Mcorr.pos = Mcorr_(:,:,index_pos);
    Mcorr.ang = Mcorr_(:,:,index_ang);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    Mcorr = [];
end

% Response matrix: orbit vs. RF frequency
fprintf('   Response matrix: orbit vs. RF frequency\n'); tic;
Mrf_ = respmrf(THERING, markers.orbit);
Mrf.pos = Mrf_(:,:,index_pos);
Mrf.ang = Mrf_(:,:,index_ang);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);

% Response matrix: orbit vs. magnet's displacements
if isfield(markers, 'disp') && ~isempty(markers.disp)
    fprintf('   Response matrix: orbit vs. magnets'' displacements\n'); tic;
    Mdisp_ = respmdisp(THERING, markers.orbit, markers.disp, cellstr(repmat(plane,length(markers.disp),1)));
    Mdisp.pos = Mdisp_(:,:,index_pos);
    Mdisp.ang = Mdisp_(:,:,index_ang);
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    Mdisp = [];
end

% Response matrix: orbit vs. ID disturbance
if isfield(markers, 'ss') && ~isempty(markers.ss)
    fprintf('   Response matrix: orbit vs. ID disturbance\n'); tic;

    nss = length(markers.ss);
    
    original_PassMethod = cell(nss,1);
    for i=1:nss
        original_PassMethod{i} = THERING{markers.ss(i)}.PassMethod;
        THERING{markers.ss(i)}.PassMethod = 'CorrectorPass';
        THERING{markers.ss(i)}.KickAngle(1) = 0;
        THERING{markers.ss(i)}.KickAngle(2) = 0;
    end
    
    Mss_ = respmcorr(THERING, markers.orbit, markers.ss, plane);
    Mss.pos = Mss_(:,:,index_pos);
    Mss.ang = Mss_(:,:,index_ang);
    
    for i=1:nss
        THERING{markers.ss(i)}.PassMethod = original_PassMethod;
    end
            
    fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);
else
    Mss = [];
end