fclog_filename = '\\stnls02.lnls.br\CorrecaoOrbita\fc\fcsend_log.txt';
ip = '10.0.5.31';
Ts = 320e-6;
npts_packet = 100;
start_from = 1;
nominal_beam_energy = 1370;

load(fullfile('test_data','temp_matrices.mat'))

respm_profiles = [];%diag(corr_steps);
bumps = (blkdiag(pinv(Mh),pinv(Mv))*([[zeros(1,16) -1 1 zeros(1,32)]; [zeros(1,16) 1 1 zeros(1,32)]; [zeros(1,41) -1 1 zeros(1,7)]; [zeros(1,41) 1 1 zeros(1,7)]])')';

if false % FIXME
    [hbpm,hsv,hcorr] = svd(Mh);
    [vbpm,vsv,vcorr] = svd(Mv);
    
    invhsv = zeros(size(hsv'));
    invhsv(1:min(size(hsv)), 1:min(size(hsv))) = diag(1./diag(hsv));
    invvsv = zeros(size(vsv'));
    invvsv(1:min(size(vsv)), 1:min(size(vsv))) = diag(1./diag(vsv));
    
    profh = invhsv'*hcorr';
    profv = invvsv'*vcorr';
    
    profh = profh(1:min(size(profh)), 1:min(size(profh)))*0.5;
    profv = profv(1:min(size(profv)), 1:min(size(profv)))*0.5;
    
    bpm_singular_vectors = blkdiag(profh, profv);
else
    bpm_singular_vectors = [];
end

manual_settings_array = {...
    %{'beam energy = 500 MeV', 500}, ...
    %{'beam energy = 1.37 GeV, IDs = open, SCW = 0.2 T', 1370}, ...
    {'AWG01 = closed, AON11 gap = 22 mm, AON11 phase = -25 mm, SCW = 3.5 T, user condition (correct tune, coupling, etc.)', 1370}, ...
    %{'SCW = 0.2 T, beam current = 30 mA', 1370}, ...
    %{'LT + booster + SR sexts = off', 1370}, ...
    };

nrepetitions = 1;
exp_descriptions_array = {...
    %{'step', 0.1}, ...
    {'ramp', 0.5}, ...
    %{'prbs', 0.1, [0 1/3]}, ...
    %{'sine', 0.1, [10 500]*2*Ts, [30 2 1]}, ...
    };

exp_descriptions = {};
ii = 1;
for i=start_from:length(manual_settings_array)
    exp_descriptions{ii} = {'manual_settings', manual_settings_array{i}{1}};
    ii=ii+1;
    for j=1:length(exp_descriptions_array)
        exp_descriptions{ii} = exp_descriptions_array{j};
        beam_energy(ii) = manual_settings_array{i}{2};
        ii=ii+1;
    end
end

for j=1:length(exp_descriptions)
    if strcmpi(exp_descriptions{j}{1}, 'manual_settings')
        if strcmpi(ip, '10.0.5.31')
            fileid = fopen(fclog_filename, 'a+');
            fprintf(fileid, ['[' fatimestr(now) '] ']);
            fprintf(fileid, ['manual settings; ' exp_descriptions{j}{2}]);
            fprintf(fileid, '\n');
            fclose(fileid);
        end
        input(sprintf('Set %s and press ENTER\n', exp_descriptions{j}{2}), 's');
    else
        expinfo = [];
        
        expinfo.excitation = exp_descriptions{j}{1};
        expinfo.amplitude = exp_descriptions{j}{2}*beam_energy(j)/nominal_beam_energy;
        
        expinfo.Ts = Ts;
        expinfo.ncols = 50;
        expinfo.duration = 100;
        expinfo.nperiods = 1;
        expinfo.pauselength = 10;
        expinfo.mode = 'corr_sum';
        expinfo.profiles = [respm_profiles; bpm_singular_vectors; bumps];
        
        if strcmpi(expinfo.excitation, 'prbs')
            expinfo.duration = 1000;
            expinfo.prbsperiod = 45;
        end
        
        if strcmpi(expinfo.excitation, 'sine')
            expinfo.sinedata = exp_descriptions{j}{4};
        end
        
        if strcmpi(expinfo.excitation, 'step') || strcmpi(expinfo.excitation, 'ramp')
            expinfo.uncorrelated = false;
            expinfo.band = [];
        else
            expinfo.uncorrelated = true;
            expinfo.band = exp_descriptions{j}{3};
        end

        stopat = (expinfo.duration+expinfo.pauselength)*(size(expinfo.profiles, 1) + 1)*nrepetitions;
        if strcmpi(ip, '10.0.5.31')
            fclog(fclog_filename, expinfo, npts_packet, stopat);
        end
        [fcdata, expout] = fcsend(ip, @fcident, expinfo, npts_packet, stopat);
    end
end
pause(10); % let data on buffer be consumed