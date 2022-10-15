function M = posang(section)
% POSANG Position and angle matrix.
% 
% Translate desired angles and positions of the bump at light source
% positions to orbit distortions at neighbor BPMs.
% From (POSX, ANGX, POSY, ANGY) to (POSX_BPM1, POSX_BPM2, POSY_BPM1, POSY_BPM2)
%
% posang(SECTION) SECTION must be 'BC', 'SA', 'SB' or 'SP', the source points.
% For 'BC', BPM1 = C2 and BPM2 = C3-1.
% For all straight sections (SA, SB and SP), BPM1 = M1 and BPM2 = M2.
%
% Sirius storage ring model SI.25.04

if strcmp(section, 'BC')
    M = [1.11371 -0.61624 0 0;
        1.25316 1.54265 0 0;
        0 0 0.90631 -0.57170;
        0 0 0.69528 1.80929];
else
    switch section
        case 'SA'
            len = 7.0358;            
        case 'SB'
            len = 6.1758;
        case 'SP'
            len = 6.1758;
    end
    
    M = [1 -len/2 0 0;
        1 len/2 0 0;
        0 0 1 -len/2;
        0 0 1 len/2];
end