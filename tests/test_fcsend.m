ip = '10.0.5.321';
Ts = 320e-6;
npts_packet = 100;
start_from = 1;

load(fullfile('test_data','temp_matrices.mat'))

if false
    [hbpm,hsv,hcorr] = svd(Mh);
    [vbpm,vsv,vcorr] = svd(Mv);
    
    invhsv = zeros(size(hsv'));
    invhsv(1:min(size(hsv)), 1:min(size(hsv))) = diag(1./diag(hsv));
    invvsv = zeros(size(vsv'));
    invvsv(1:min(size(vsv)), 1:min(size(vsv))) = diag(1./diag(vsv));
    
    profh = invhsv'*hcorr';
    profv = invvsv'*vcorr';
    
    profh = profh(1:min(size(profh)), 1:min(size(profh)))*0.5;
    profv = profv(1:min(size(profv)), 1:min(size(profv)))*0.05;
    
    bpm_singular_vectors = blkdiag(profh, profv);
else
    bpm_singular_vectors = [];
end

manual_settings_array = {...
    'beam energy = 500 MeV', ...
    'beam energy = 900 MeV', ...
    'beam energy = 1.37 GeV', ...
    'IDs = open, SCW = 0.2 T', ...
    'AWG01 = closed, AON11 = closed', ...
    'SCW = 2 T', ...
    'SCW = 3.8 T, user condition (correct tune, coupling, etc.)', ...
    'SCW = 0.2 T, beam current = 30 mA', ...
    'LT + booster + SR sexts = off', ...
};

nrepetitions = 1;
exp_descriptions_array = {...
    {'step', 1}, ...
    {'prbs', 1, [0 0.25]}, ...
    {'prbs', 1, [0 0.5]}, ...
    {'prbs', 1, [0 1]}, ...
    {'sine', 1, [10 500]*2*Ts, [30 2 1]}, ...
    {'rgs',  1, [0 1]}, ...
    {'step', 0.5}, ...
    {'prbs', 0.5, [0 0.25]}, ...
    {'prbs', 0.5, [0 0.5]}, ...
    {'prbs', 0.5, [0 1]}, ...
    {'sine', 0.5, [10 500]*2*Ts, [30 2 1]}, ...
    {'rgs',  0.5, [0 1]}, ...
    {'step', 0.25}, ...
    {'prbs', 0.25, [0 0.25]}, ...
    {'prbs', 0.25, [0 0.5]}, ...
    {'prbs', 0.25, [0 1]}, ...
    {'sine', 0.25, [10 500]*2*Ts, [30 2 1]}, ...
    {'rgs',  0.25, [0 1]}, ...
};

exp_descriptions = {};
ii = 1;
for i=start_from:length(manual_settings_array)
    exp_descriptions{ii} = {'manual_settings', manual_settings_array{i}};
    ii=ii+1;
    for j=1:length(exp_descriptions_array)
        exp_descriptions{ii} = exp_descriptions_array{j};
        ii=ii+1;
    end
end

for j=1:length(exp_descriptions)
    if strcmpi(exp_descriptions{j}{1}, 'manual_settings')
        if strcmpi(ip, '10.0.5.31')
            fileid = fopen(filename, 'a+');
            fprintf(fileid, ['[' fatimestr(now) '] ']);
            fprintf(fileid, ['manual settings; ' exp_descriptions{j,2}]);
            fprintf(fileid, '\n');
            fclose(fileid);
        end
        input(sprintf('Set %s and press ENTER\n', exp_descriptions{j}{2}), 's');
    else
        expinfo = [];
        
        expinfo.excitation = exp_descriptions{j}{1};
        expinfo.amplitude = exp_descriptions{j}{2};
        
        if strcmpi(expinfo.excitation, 'sine')
            expinfo.sinedata = exp_descriptions{j}{4};
        end
        
        if strcmpi(expinfo.excitation, 'step')
            expinfo.uncorrelated = false;
        else
            expinfo.uncorrelated = true;
            expinfo.band = exp_descriptions{j}{3};
        end
        
        expinfo.Ts = Ts;
        expinfo.ncols = 50;
        expinfo.duration = 100;
        expinfo.pauselength = 10;
        expinfo.mode = 'corr_sum';
        expinfo.profiles = [diag(corr_steps); bpm_singular_vectors];

        stopat = (expinfo.duration+expinfo.pauselength)*(size(expinfo.profiles, 1) + 1)*nrepetitions;
        if strcmpi(ip, '10.0.5.31')
            fclog('\\stnls02.lnls.br\CorrecaoOrbita\fc\fcsend_log.txt', expinfo, npts_packet, stopat);
        end
        [fcdata, expout] = fcsend(ip, @fcident, expinfo, npts_packet, stopat);
    end
end