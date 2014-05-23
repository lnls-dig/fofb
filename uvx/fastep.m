function [dbpm, dcorr] = fastep(fadata, ps_index, threshold)

index1 = find(fadata.ps_readings(:,ps_index) < threshold);
index2 = find(fadata.ps_readings(:,ps_index) >= threshold);

bpm1 = mean(fadata.bpm_readings(index1,:));
bpm2 = mean(fadata.bpm_readings(index2,:));
dbpm = bpm2-bpm1;

corr1 = mean(fadata.ps_readings(index1,:));
corr2 = mean(fadata.ps_readings(index2,:));
dcorr = corr2-corr1;