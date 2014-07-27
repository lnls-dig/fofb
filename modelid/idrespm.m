function M = idrespm(bpm_iddata, npts_level)

corr_setpoints = bpm_iddata{1}.InputData;
bpm_readings = bpm_iddata{1}.OutputData;
npts = size(bpm_readings, 1);

orb_before_step = zeros(size(bpm_readings,2),size(corr_setpoints,2));
orb_after_step = zeros(size(bpm_readings,2),size(corr_setpoints,2));
corr_steps = zeros(1,size(corr_setpoints,2));

for i=1:length(bpm_iddata)-1
    corr_setpoints = bpm_iddata{i+1}.InputData;
    bpm_readings = bpm_iddata{i+1}.OutputData;
    npts = size(bpm_readings, 1);
    
    orb_before_step(:,i) = levelmeanstd(bpm_readings, npts_level, npts/2);
    orb_after_step(:,i) = levelmeanstd(bpm_readings, npts_level, npts);

    corr_steps(1,i) = corr_setpoints(end) - corr_setpoints(end/2);
end

M = double((orb_after_step-orb_before_step)./repmat(corr_steps, size(orb_before_step, 1), 1));


function [level_mean, level_std] = levelmeanstd(data, npts_level, idx)

level = data(idx+(-npts_level+1:0),:);
level_mean = mean(level);
level_std = std(level);