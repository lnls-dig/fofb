function r = respmuvx

global THERING;

% Orbit and light source markers
bpm    = sort(findcells(THERING, 'FamName', 'BPM'))';
rad4g  = findcells(THERING, 'FamName', 'RAD4G')';
rad15g = findcells(THERING, 'FamName', 'RAD15G')';
ss     = findcells(THERING, 'FamName', 'LCENTER')';
id = sort(ss([1 5:end])); % Remove injection (2) and cavity (3) long straight sections; remove long straight section without ID (4)
light_source = sort([rad4g; rad15g; id]);

% Actuator and disturbance markers
hcm = findcells(THERING, 'FamName', 'HCM')';
vcm = findcells(THERING, 'FamName', 'VCM')';

nsegsquad = 2;
quad = ([ ...
    findcells(THERING, 'FamName', 'A2QF01') ...
    findcells(THERING, 'FamName', 'A2QF03') ...
    findcells(THERING, 'FamName', 'A2QF05') ...
    findcells(THERING, 'FamName', 'A2QF07') ...
    findcells(THERING, 'FamName', 'A2QF09') ...
    findcells(THERING, 'FamName', 'A2QF11') ...
    findcells(THERING, 'FamName', 'A2QD01') ...
    findcells(THERING, 'FamName', 'A2QD03') ...
    findcells(THERING, 'FamName', 'A2QD05') ...
    findcells(THERING, 'FamName', 'A2QD07') ...
    findcells(THERING, 'FamName', 'A2QD09') ...
    findcells(THERING, 'FamName', 'A2QD11') ...
    findcells(THERING, 'FamName', 'A6QF01') ...
    findcells(THERING, 'FamName', 'A6QF02') ...
])';
[~, reorder_quad_markers] = sort(quad(1:nsegsquad:end));
quad = reshape(quad, nsegsquad, length(quad)/nsegsquad)';
quad = quad(reorder_quad_markers, :);

% % Orbit and light source markers
% bpm_names    = {'01B' '02A' '02B' '03A' '03B' '03C' '04A' '04B' '05A' '05B' '06A' '06B' '07A' '07B' '08A' '08B' '09A' '09B' '10A' '10B' '11A' '11B' '12A' '12B' '01A'};
% rad4g_names  = {'D01A' 'D02A (SAXS2)' 'D03A (IR)' 'D04A (SXS)' 'D05A (TGM)' 'D06A (DXAS)' 'D07A' 'D08A (SGM)' 'D09A' 'D10A (XRD2)' 'D11A' 'D12A (XRD1)' };
% rad15g_names = {'D01B (SAXS1)' 'D02B (DFX)' 'D03B (MX1)' 'D04B (XAFS1)' 'D05B (DFE)' 'D06B (IMX1)' 'D07B' 'D08B (XAFS2)' 'D09B (XRF)' 'D10B (XPD)' 'D11B' 'D12B'};
% ss_names     = {'AWG01', 'INJECTION', 'RF', '(empty SS)', 'AWG09', 'AON11'};
% id_names     = {'W01B (MX2)', 'W09A (XDS)', 'U11A (PGM1)'};
% source_names = [rad4g_names rad15g_names id_names];
% source_names = source_names(reorder_source_markers);
% orbit_names  = [bpm_names source_names];
% orbit_names  = orbit_names(reorder_orbit_markers);
% % Actuator and disturbance names
% hcm_names    = {'ACH01B' 'ACH02' 'ACH03A' 'ACH03B' 'ACH04' 'ACH05A' 'ACH05B' 'ACH06' 'ACH07A' 'ACH07B' 'ACH08' 'ACH09A' 'ACH09B' 'ACH10' 'ACH11A' 'ACH11B' 'ACH12' 'ACH01A'};
% vcm_names    = {'ACV01B' 'ALV02A' 'ALV02B' 'ACV03A' 'ACV03B' 'ALV04A' 'ALV04B' 'ACV05A' 'ACV05B' 'ALV06A' 'ALV06B' 'ACV07A' 'ACV07B' 'ALV08A' 'ALV08B' 'ACV09A' 'ACV09B' 'ALV10A' 'ALV10B' 'ACV11A' 'ACV11B' 'ALV12A' 'ALV12B' 'ACV01A'};
% quad_names   = {'A2QF01B' 'A2QF01A' 'A2QF03A' 'A2QF03B' 'A2QF05A' 'A2QF05B' 'A2QF07A' 'A2QF07B' 'A2QF09A' 'A2QF09B' 'A2QF11A' 'A2QF11B' ...
%                 'A2QD01B' 'A2QD01A' 'A2QD03A' 'A2QD03B' 'A2QD05A' 'A2QD05B' 'A2QD07A' 'A2QD07B' 'A2QD09A' 'A2QD09B' 'A2QD11A' 'A2QD11B' ...
%                 'A6QF02A' 'A6QF04B' 'A6QF06A' 'A6QF08B' 'A6QF10A' 'A6QF12B' ...
%                 'A6QF02B' 'A6QF04A' 'A6QF06B' 'A6QF08A' 'A6QF10B' 'A6QF12A'};
% quad_names   = quad_names(reorder_quad_markers);
% rf_names     = {'RF A', 'RF B'};

% Calculate response matrices
[r.bpm.fofb.Mx, r.bpm.fofb.Mx_] = respmcorr(bpm, hcm, 'x');
[r.bpm.fofb.My, r.bpm.fofb.My_] = respmcorr(bpm, vcm, 'y');
[r.bpm.dist.id.Mx, r.bpm.dist.id.Mx_] = respmdisp(bpm, id, 'x');
[r.bpm.dist.id.My, r.bpm.dist.id.My_] = respmdisp(bpm, id, 'y');
[r.bpm.dist.misalign.quad.Mx, r.bpm.dist.misalign.quad.Mx_] = respmdisp(bpm, quad , 'x');
[r.bpm.dist.misalign.quad.My, r.bpm.dist.misalign.quad.My_] = respmdisp(bpm, quad , 'y');
[r.bpm.dist.rf.M, r.bpm.dist.rf.M_] = respmrf(bpm);

[r.light_source.fofb.Mx, r.light_source.fofb.Mx_] = respmcorr(light_source, hcm, 'x');
[r.light_source.fofb.My, r.light_source.fofb.My_] = respmcorr(light_source, vcm, 'y');
[r.light_source.dist.id.Mx, r.light_source.dist.id.Mx_] = respmdisp(light_source, id, 'x');
[r.light_source.dist.id.My, r.light_source.dist.id.My_] = respmdisp(light_source, id, 'y');
[r.light_source.dist.misalign.quad.Mx, r.light_source.dist.misalign.quad.Mx_] = respmdisp(light_source, quad , 'x');
[r.light_source.dist.misalign.quad.My, r.light_source.dist.misalign.quad.My_] = respmdisp(light_source, quad , 'y');
[r.light_source.dist.rf.M, r.light_source.dist.rf.M_] = respmrf(light_source);