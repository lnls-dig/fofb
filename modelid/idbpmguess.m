% ach = 1:18;
% acv = sort([19:4:42 22:4:42]);
% alv = setdiff(19:42,acv);

prbsperiod = 2047;

order_corr = [4 3 0];
order_bpm = [2 2 1];

bpmhvsymmetry = 0;
bpmi = 1:50;
corri = 1:42;

bpmungain = 1;

datae = {};
datav = {};
for i=1:length(corri)
    datae{i} = corr_iddata_id{corri(i)+1};
    datav{i} = corr_iddata_val{corri(i)+1};
end

fixedsys = tf([0.18 0.64 0.18], [1 0 0],  -1, 'InputDelay', 1);
[acorr, bcorr, fitcorr] = idsiso(datae, datav, order_corr, 1, fixedsys);

syscorrps = idab2sys(acorr,bcorr);

% --

bpmi = bpmi(:);

M_ = M(bpmi, corri);
[~, corribpmi] = max(abs(M_), [], 2);
gains = M_(sub2ind(size(M_),(1:size(M_,1))',corribpmi));

datae = {};
datav = {};
for i=1:length(bpmi)
    datae{i} = bpm_iddata_id{corri(corribpmi(i))+1}(:,bpmi(i));
    datav{i} = bpm_iddata_val{corri(corribpmi(i))+1}(:,bpmi(i));
    
    datae{i}.OutputData = datae{i}.OutputData/gains(i);
    datav{i}.OutputData = datav{i}.OutputData/gains(i);
end
[abpm, bbpm, fitbpm] = idsiso(datae, datav, order_bpm, syscorrps(corribpmi), fixedsys);

sysbpm = idab2sys(abpm,bbpm);

if bpmhvsymmetry
    % Take only vertical BPM responses
    abpm = abpm(end/2+1:end, :);
    bbpm = bbpm(end/2+1:end, :);
end

% --

[datae, datav] = idpreprocessdata(bpm_iddata_id{1}, bpm_iddata_val{1}, corri, bpmi, 1, fixedsys, prbsperiod);

datae.InputData = datae.InputData*1e4;
datae.OutputData = datae.OutputData*1e4;
datav.InputData = datav.InputData*1e4;
datav.OutputData = datav.OutputData*1e4;

[param, fitfofb] = idgbbpm(datae, datav, order_bpm, order_corr, abpm, bbpm, acorr, bcorr, M_, bpmungain, bpmhvsymmetry);


pbpmoffset = 0;
if ~bpmungain
    pcorroffset = size(abpm,1)*(sum(order_bpm(1:2)));
else
    pcorroffset = size(abpm,1)*(sum(order_bpm(1:2))-1);
end

parambpm = param(pbpmoffset+1:pcorroffset);
paramcorr = param(pcorroffset+1:end);

[abpm_idgb,bbpm_idgb] = idparam2ab(parambpm, size(abpm,1), order_bpm, bpmungain);
[acorr_idgb,bcorr_idgb] = idparam2ab(paramcorr, size(acorr,1), order_corr);