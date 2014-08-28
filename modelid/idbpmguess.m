% ach = 1:18;
% acv = sort([19:4:42 22:4:42]);
% alv = setdiff(19:42,acv);

prbsperiod = 2047;

if false
    load(fullfile('~', 'dados_ident','master','sysident_data_2014.07.27','processed','beam energy = 1.37 GeV, IDs = closed','AON11 gap = 22 mm, SCW = 3.5 T, user condition (correct tune, coupling, etc.)','step','2014.07.27_17.23.42.mat'))
    M=idrespm(bpm_iddata, 6000);
    
    load(fullfile('~', 'dados_ident','master','sysident_data_2014.07.27','processed','beam energy = 1.37 GeV, noise','2014.07.27.mat'));
    
    load(fullfile('~', 'dados_ident','master','sysident_data_2014.07.27','processed','beam energy = 1.37 GeV, IDs = closed','AON11 gap = 22 mm, SCW = 3.5 T, user condition (correct tune, coupling, etc.)','prbs','2014.07.27_17.28.30.mat'))
    [bpm_iddata_id, bpm_iddata_val] = idprbsmean(bpm_iddata, prbsperiod, [], false);
    [corr_iddata_id, corr_iddata_val] = idprbsmean(corr_iddata, prbsperiod, [], false);
end

%idmethod_corrps = @n4sid;
idmethod_corrps = @arx;
order_corrps = [4 4 0];

% idmethod_bpm = @n4sid;
idmethod_bpm = @arx;
order_bpm = [2 2 1];

bpmi = 1:50;
corri = 1:42;

bpmungain = 1;
corrpsungain = 0;
bpmhvsymmetry = 0;

cicfilter_netdelay = tf([0.18 0.64 0.18], [1 0 0],  -1, 'InputDelay', 1);

Ts = bpm_iddata_id{1}.Ts;


nbpm = length(bpmi);
ncorr = length(corri);

nnoise = 100;

% for i=1:ncorr
%     tic
%     m = ar(corr_iddata_noise.OutputData(:,i), nnoise);
%     syscorrnoise{i} = m.a;
%     toc
% end
% 
% for i=1:nbpm
%     tic
%     m = ar(bpm_iddata_noise.OutputData(:,i), nnoise);
%     sysbpmnoise{i} = m.a;
%     toc
% end


% Black box identification of orbit corrector power supplies models
datae = {};
datav = {};
for i=1:ncorr
    datae{i} = idfilt(detrend(corr_iddata_id{corri(i)+1},0), {syscorrnoise{i},[1 zeros(1,length(syscorrnoise{i}))]}, 'causal');
    datav{i} = idfilt(detrend(corr_iddata_val{corri(i)+1},0), {syscorrnoise{i},[1 zeros(1,length(syscorrnoise{i}))]}, 'causal');
end

[acorrps, bcorrps, fitcorrps, rcorrps] = idsiso(datae, datav, idmethod_corrps, order_corrps, 1, cicfilter_netdelay);
syscorrps = idab2sys(acorrps, bcorrps, Ts);

% Black box identification of BPM electronics + beam dynamics + vacuum chamber + orbit corrector magnetic core model
% Use experiments where only one orbit corrector was excited. For each BPM, take the corrector that produces the largest distortion at the location of that BPM.
M_ = M(bpmi, corri);
[~, corribpmi] = max(abs(M_), [], 2);
gains = M_(sub2ind(size(M_),(1:size(M_,1))',corribpmi));

datae = {};
datav = {};
for i=1:nbpm
    datae{i} = idfilt(detrend(bpm_iddata_id{bpmi(i)+1},0), {sysbpmnoise{i},[1 zeros(1,length(sysbpmnoise{i}))]}, 'causal');
    datav{i} = idfilt(detrend(bpm_iddata_val{bpmi(i)+1},0), {sysbpmnoise{i},[1 zeros(1,length(sysbpmnoise{i}))]}, 'causal');
    
    datae{i}.OutputData = datae{i}.OutputData/gains(i);
    datav{i}.OutputData = datav{i}.OutputData/gains(i);
end

[abpm, bbpm, fitbpm, rbpm] = idsiso(datae, datav, idmethod_bpm, order_bpm, syscorrps(corribpmi), cicfilter_netdelay);
sysbpm = idab2sys(abpm, bbpm, Ts);

% Grey box identification of the complete FOFB model using individual orbit corrector, BPM and static beam response matrix
% as starting point and static beam response matrix obtained from step experiments
[datae, datav] = idpreprocessdata(bpm_iddata_id{1}, bpm_iddata_val{1}, corri, bpmi, 1, cicfilter_netdelay);

% Scaling factor applied to input-output data in order to avoid reaching floating point resolution limits on cost function calculation 
f = 1e2;
datae.InputData = datae.InputData*f;
datae.OutputData = datae.OutputData*f;
datav.InputData = datav.InputData*f;
datav.OutputData = datav.OutputData*f;

if bpmhvsymmetry
    % Take only vertical BPM responses
    abpm = abpm(end/2+1:end, :);
    bbpm = bbpm(end/2+1:end, :);
end

num = num2cell(zeros(ncorr));
den = num2cell(ones(ncorr));
for i=1:ncorr, num{i,i} = bcorrps{i}; den{i,i} = acorrps{i}; end
syscorrpsmimo = tf(num,den,320e-6);

num = num2cell(zeros(nbpm));
den = num2cell(ones(nbpm));
for i=1:nbpm, num{i,i} = bbpm{i}; den{i,i} = abpm{i}; end
sysbpmmimo = tf(num,den,320e-6);

[param, fitfofb] = idgbbpm(datae, datav, order_bpm, order_corrps, abpm, bbpm, acorrps, bcorrps, M_, bpmungain, corrpsungain, bpmhvsymmetry);

% Presentation of results
pbpmoffset = 0;
if ~bpmungain
    pcorroffset = size(abpm,1)*(sum(order_bpm(1:2)));
else
    pcorroffset = size(abpm,1)*(sum(order_bpm(1:2))-1);
end

param_bpm = param(pbpmoffset+1:pcorroffset);
param_corrps = param(pcorroffset+1:end);

[abpm_idgb,bbpm_idgb] = idparam2ab(param_bpm, size(abpm,1), order_bpm, bpmungain);
[acorr_idgb,bcorr_idgb] = idparam2ab(param_corrps, size(acorrps,1), order_corrps, corrpsungain);

sysbpm2 = idab2sys(abpm_idgb, bbpm_idgb, Ts);
syscorrps2 = idab2sys(acorr_idgb, bcorr_idgb, Ts);