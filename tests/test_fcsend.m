load(fullfile('test_data','temp_matrices.mat'))
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

npts_packet = 100;

expinfo.Ts = 320e-6;
expinfo.ncols = 50;
expinfo.excitation = 'sine';
expinfo.amplitude = 0.2;
expinfo.duration = 10;
expinfo.pauselength = 10;
expinfo.mode = 'corr_sum';
expinfo.profiles = [diag(corr_steps); bpm_singular_vectors];

if strcmpi(expinfo.excitation, 'step')
    expinfo.uncorrelated = false;
else
    expinfo.uncorrelated = true;
end

if strcmpi(expinfo.excitation, 'sine')
    expinfo.sinedata = [10 2 1];
    expinfo.band = [10 500]*expinfo.Ts*2;
elseif strcmpi(expinfo.excitation, 'rgs') || strcmpi(expinfo.excitation, 'rbs') || strcmpi(expinfo.excitation, 'prbs')
    expinfo.band = [0 0.5];
end

stopat = (expinfo.duration+expinfo.pauselength)*(size(expinfo.profiles, 1) + 1);
fclog('\\stnls02.lnls.br\CorrecaoOrbita\fc\fcsend_log.txt', expinfo, npts_packet, stopat);
[fcdata, expout] = fcsend('10.0.5.31', @fcident, expinfo, npts_packet, stopat);