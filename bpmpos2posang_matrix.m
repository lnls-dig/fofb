function [posx, posy, angx, angy] = bpmpos2posang_matrix

Msa = inv(posang('SA'));
Msb = inv(posang('SB'));
Msp = inv(posang('SP'));
Mbc = inv(posang('BC'));

% Inv. posang translates:
% [bpm1_posx, bpm2_posx, bpm1_posy, bpm2_posy] into [posx,angx,posy,angy]

Z = zeros(1,2);

for i=1:2 %[Horiozntal, Vertical] Position or Angle at light source
    for j=1:2 % [Horizontal, Vertical] BPM
        for k=1:2 %[Position, Angle] at light source
            M_ = blkdiag(...
                Msa(2*(i-1)+k, [2*j-1 2*j]), Z, Mbc(2*(i-1)+k, [2*j-1 2*j]), Z, ...
                Msb(2*(i-1)+k, [2*j-1 2*j]), Z, Mbc(2*(i-1)+k, [2*j-1 2*j]), Z, ...
                Msp(2*(i-1)+k, [2*j-1 2*j]), Z, Mbc(2*(i-1)+k, [2*j-1 2*j]), Z, ...
                Msb(2*(i-1)+k, [2*j-1 2*j]), Z, Mbc(2*(i-1)+k, [2*j-1 2*j]), Z);
            M_ = blkdiag(M_,M_,M_,M_,M_);
            % Take care of first M1 BPM of the ring, which is the one of
            % higher index according to the ring convention (last BPM to
            % "see" an injected beam)
            M_ = [M_(:,2:end) M_(:,1)];
            M{i,j,k} = M_;
        end
    end
end

posx.bpmx = M{1,1,1};
posx.bpmy = M{1,2,1};
posy.bpmx = M{2,1,1};
posy.bpmy = M{2,2,1};

angx.bpmx = M{1,1,2};
angx.bpmy = M{1,2,2};
angy.bpmx = M{2,1,2};
angy.bpmy = M{2,2,2};