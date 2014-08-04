Ts = 320e-6;
prbsperiod=2047;
ncorrmod=4;
nk = 1;
nu = 1;
mcorr = idgrey('idcorrpscicmodel',zeros(nu*(2*ncorrmod+1),1),'d',[nk,ncorrmod,nu],'InitialState','Zero');
mcorr.Ts = Ts;
%'DisturbanceModel','None'

syscorrarray = {};
for corr=1:42
    corr_iddata_id{corr+1}.Ts = Ts;
    m = tf(pem(detrend(corr_iddata_id{corr+1},0), mcorr));
    syscorrarray{corr} = m(1,1);
end

bpms = 1:25;
corrs = 1:18;

nbpm=length(bpms);
ncorr=length(corrs);
nabpm=3;
nbbpm=2;
nacorr=0;
nbcorr=0;
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
param = [param; zeros(nbpm*(nabpm+nbbpm),1)];
%param = [param; randn(nbpm,1)];
% param = [param; ...
%      repmat([ -0.0832876040803259; 0.655750344618449; -0.13946289888987], nbpm, 1); ...
%      repmat([ 0.0349942423480944; 0.292435068837765; 0.205868964777135], nbpm, 1)];
%load parameters;
%param = paramhv;
if ~bpmungain
    param = [param; repmat(0.0193741179554821, nbpm, 1)];
end

param = [param; zeros(ncorr*(nacorr+nbcorr),1)];
% param = [param; ones(ncorr,1)];
if ~corrungain
    param = [param; ones(ncorr,1)];
end

auxvars = [nbpm,ncorr,nabpm,nbbpm,nacorr,nbcorr,bpmungain,corrungain,M_(:)'];
fofbmodel = idgrey('idfofbmodel',param,'d', auxvars,'InitialState','Zero');
fofbmodel.InputDelay = repmat(nk, ncorr, 1);
fofbmodel.Ts = Ts;
fofbmodel.Algorithm.Focus = 'Simulation';
fofbmodel.Algorithm.Display = 'On';

datae.Ts = Ts;
fofbmodele = pem(detrend(datae,0), fofbmodel);

[~,fit] = compare(fofbmodele, detrend(datav,0), Inf);
fit = squeeze(fit);
hold all; plot(fit)
