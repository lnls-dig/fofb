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

bpms = 1:25;
corrs = 1:18;

nbpm=length(bpms);
ncorr=length(corrs);
nxbpm=2;
nxcorr=0;
nk=1;

M_ = M(bpms,corrs);

ue=[];
ye=[];
uv=[];
yv=[];

i=1;
for corr = corrs
    u = bpm_iddata_id{1}.InputData(:,corr);
    um = lsim(syscorrarray{corr}, [u; u]);
    ue(:,i) = um(prbsperiod+1:2*prbsperiod);
    uv(:,i) = um(prbsperiod+1:2*prbsperiod);
    i=i+1;
end

%param = 

param = zeros(nbpm*(2*nxbpm+1) + ncorr*(2*nxcorr+1),1);
%load parameters; param = repmat(param, 1, 2);

param = [repmat([-0.0399 0.5288], 1, nbpm)  repmat([0.6777 0.1918], 1, nbpm) repmat(1, 1, nbpm)];

fofbmodel = idgrey('idfofbmodel',param,'d',[nk,nbpm,ncorr,nxbpm,nxcorr,M_(:)'],'InitialState','Zero');
fofbmodel.Ts = Ts;
fofbmodel.Algorithm.Focus = 'Simulation';
fofbmodel.Algorithm.Display = 'On';

datae = iddata(bpm_iddata_id{1}.OutputData(:,bpms),ue,Ts);
datav = iddata(bpm_iddata_id{1}.OutputData(:,bpms),uv,Ts);

datae.Ts = Ts;
m4 = pem(detrend(datae,0), fofbmodel);

[~,fit] = compare(m4,detrend(datav,0));
fit = squeeze(fit);
hold all; plot(fit)
