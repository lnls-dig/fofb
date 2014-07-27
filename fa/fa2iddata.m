function [corr_iddata, bpm_iddata] = fa2iddata(fadata)

corr_iddata = iddata( ...
    'InputData', double(fadata.corr_setpoints), ...
    'OutputData', double(fadata.corr_readings), ...
    'Ts', repmat(double(fadata.period)*1e-6, size(fadata.corr_setpoints, 2), 1), ...
    'InputName', fadata.corr_setpoints_names, ...
    'OutputName', fadata.corr_names, ...
    'InputUnit', repmat({'A'}, size(fadata.corr_readings, 2), 1), ...
    'OutputUnit', repmat({'A'}, size(fadata.corr_setpoints, 2), 1), ...
    'TimeUnit', 's' ...
    );

bpm_iddata = iddata( ...
    'InputData', double(fadata.corr_setpoints), ...
    'OutputData', double(fadata.bpm_readings), ...
    'Ts', repmat(double(fadata.period)*1e-6, size(fadata.corr_setpoints, 2), 1), ...
    'InputName', fadata.corr_setpoints_names, ...
    'OutputName', fadata.bpm_names, ...
    'InputUnit', repmat({'A'}, size(fadata.corr_readings, 2), 1), ...
    'OutputUnit', repmat({'mm'}, size(fadata.bpm_readings, 2), 1), ...
    'TimeUnit', 's' ...
    );