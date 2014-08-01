% Ts = 320e-6;
% prbsperiod=2047;
% ncorrmod=4;
% nk = 1;
% nu = 1;
% mcorr = idgrey('idcorrpscicmodel',zeros(nu*(2*ncorrmod+1),1),'d',[nk,ncorrmod,nu],'InitialState','Zero');
% mcorr.Ts = Ts;
% %'DisturbanceModel','None'
% 
% syscorrarray = {};
% for corr=1:42
%     corr_iddata_id{corr+1}.Ts = Ts;
%     m3 = tf(pem(detrend(corr_iddata_id{corr+1},0), mcorr));
%     syscorrarray{corr} = m3(1,1);
% end

bpms = 1:50;
corrs = 1:42;

nbpm=length(bpms);
ncorr=length(corrs);
nxbpm=3;
nxcorr=0;
nk=0;
bpmungain = true;
corrungain = true;

M_ = M(bpms,corrs);

ue=[];
ye=[];
uv=[];
yv=[];

nexp=1;
i=1;
for corr = corrs
    if nexp == 1
        u = bpm_iddata_id{nexp}.InputData(:,corr);
    else
        if corr == nexp-1
            u = bpm_iddata_id{nexp}.InputData(:,1);
        else
            u = zeros(prbsperiod,1);
        end
    end
    um = lsim(syscorrarray{corr}, [u; u]);
    ue(:,i) = um(prbsperiod+1:2*prbsperiod);
    uv(:,i) = um(prbsperiod+1:2*prbsperiod);
    i=i+1;
end

datae = iddata(bpm_iddata_id{nexp}.OutputData(:,bpms),ue,Ts);
datav = iddata(bpm_iddata_val{nexp}.OutputData(:,bpms),uv,Ts);

param = [];
%param = [param; zeros(nbpm*(2*nxbpm),1); ones(nbpm,1)];
param = [param; ...
    repmat([0.153478424046198; 0.492048732715991; -0.361653828353039], nbpm, 1); ...
    repmat([0.318634552964485; 0.279272511234536; 0.0354099666806719], nbpm, 1)];

if ~bpmungain
    param = [param; repmat(0.0326171204852166, nbpm, 1)];
end

param = [param; zeros(ncorr*(2*nxcorr),1)];
if ~corrungain
    param = [param; ones(ncorr,1)];
end

auxvars = [nbpm,ncorr,nxbpm,nxcorr,bpmungain,corrungain,M_(:)'];
fofbmodel = idgrey('idfofbmodel',param,'d', auxvars,'InitialState','Zero');
fofbmodel.InputDelay = repmat(nk, ncorr, 1);
fofbmodel.Ts = Ts;
fofbmodel.Algorithm.Focus = 'Simulation';
fofbmodel.Algorithm.Display = 'On';

datae.Ts = Ts;
m4 = pem(detrend(datae,0), fofbmodel);

[~,fit] = compare(m4,detrend(datav,0), Inf);
fit = squeeze(fit);
hold all; plot(fit)
