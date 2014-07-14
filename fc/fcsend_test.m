expinfo.Ts = 320e-6;
expinfo.ncols = 50;

expinfo.excitation = 'step';
expinfo.amplitude = 0.1;
expinfo.duration = 100;
expinfo.mode = 'corr_sum';
expinfo.profiles = eye(42);
expinfo.uncorrelated = false;

if strcmpi(expinfo.excitation, 'sine')
    expinfo.sinedata = [10 2 1];
    expinfo.band = [10 500]*expinfo.Ts*2;
elseif strcmpi(expinfo.excitation, 'rgs') || ...
       strcmpi(expinfo.excitation, 'prbs')
    expinfo.band = [0 1];
end