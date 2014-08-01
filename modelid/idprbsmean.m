function [sysid_data_array_id, sysid_data_array_val] = idprbsmean(sysid_data_array, period)

factor = 0.5;

for i=1:length(sysid_data_array)
    sysid_data = sysid_data_array{i};
    npts = size(sysid_data, 1);
    nperiods = floor(npts/period);
    
    nperiods1 = floor(factor*(nperiods-1));
    nperiods2 = nperiods-1-nperiods1;
    
    sysid_data_id = sysid_data(period+(1:period*nperiods1));
    sysid_data_val = sysid_data(period*nperiods1+(1:period*nperiods2));
    
    for j=1:size(sysid_data_id.OutputData,2)
       sysid_data_id.OutputData(1:period,j) = mean(reshape(sysid_data_id.OutputData(:,j), period, nperiods1), 2);
       sysid_data_val.OutputData(1:period,j) = mean(reshape(sysid_data_val.OutputData(:,j), period, nperiods2), 2);
    end
    
    sysid_data_id = sysid_data_id(1:period);
    sysid_data_val = sysid_data_val(1:period);
    
    sysid_data_id.Tstart = 0;
    sysid_data_val.Tstart = size(sysid_data_id,1)*sysid_data_id.Ts;
    
    sysid_data_array_id{i} = sysid_data_id;
    sysid_data_array_val{i} = sysid_data_val;
end