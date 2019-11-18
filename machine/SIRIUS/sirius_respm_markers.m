sirius_families = sirius_si_family_data(THERING);

% BPMs, IDs and source points
bpm = sirius_families.BPM.ATIndex;
mc  = findcells(THERING, 'FamName', 'mc')';
mia = findcells(THERING, 'FamName', 'mia')';
mib = findcells(THERING, 'FamName', 'mib')';
mip = findcells(THERING, 'FamName', 'mip')';
id = sort([mia; mib; mip]);
source =  sort([mc; id]);

% Magnets
ch = sirius_families.CH.ATIndex;
cv = sirius_families.CV.ATIndex;
fc = sirius_families.FC.ATIndex;
quad = sirius_families.QN.ATIndex;
bc = sirius_families.BC.ATIndex;
b1 = sirius_families.B1.ATIndex;
b2 = sirius_families.B2.ATIndex;

% RF
rf = findcells(THERING, 'FamName', 'SRFCav')';