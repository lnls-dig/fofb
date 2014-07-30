Ts = 320e-6;
prbsperiod=2047;
ncorrmod=4;
nk = 1;
nu = 1;
mcorr = idgrey('idcorrcicmodel',zeros(nu*(2*ncorrmod+1),1),'d',[nk,ncorrmod,nu],'InitialState','Zero');
mcorr.Ts = Ts;
%'DisturbanceModel','None'

syscorrarray = {};
for corr=1:42
    corr_iddata_id{corr+1}.Ts = Ts;
    m3 = tf(pem(detrend(corr_iddata_id{corr+1},0), mcorr));
    syscorrarray{corr} = m3(1,1);
end

sysbpm_array = {};
fit = {};

nbpm=25;
ncorr=18;
nbpmmod=3;
nk=1;

for bpm=1:nbpm
    
    ue=[];
    ye=[];
    uv=[];
    yv=[];
    
    for corr=1:ncorr
        u = bpm_iddata_id{corr+1}.InputData;
        um = lsim(syscorrarray{corr}, [u; u]);
        
        ue(:,corr) = um(prbsperiod+1:2*prbsperiod);
        ye(:,corr) = bpm_iddata_id{corr+1}.OutputData(:,bpm);
        uv(:,corr) = um(prbsperiod+1:2*prbsperiod);
        yv(:,corr) = bpm_iddata_val{corr+1}.OutputData(:,bpm);
    end
    
    std_factor = std(ye)./sort(std(ue));
    selected_corrs = find(std_factor/max(std_factor) > 0.5);

    mbpm = idgrey('idbpmmodel',zeros(2*nbpmmod+1+length(selected_corrs),1),'d',[nk,nbpmmod,length(selected_corrs)],'InitialState','Zero');
    mbpm.Ts = Ts;
    
    ue = ue(:, selected_corrs);
    ye = ye(:, selected_corrs);
    uv = uv(:, selected_corrs);
    yv = yv(:, selected_corrs);
    
    datae = iddata(ye,ue,Ts);
    datav = iddata(yv,uv,Ts);
    
    datae.Ts = Ts;
    m4 = pem(detrend(datae,0), mbpm);
    
    sysbpm = ss(m4);
    sysbpm = ss(sysbpm.A(1:order,1:order), sysbpm.B(1:order,1),sysbpm.C(1,1:order),sysbpm.D(1,1),320e-6);
    sysbpm_array{bpm} = tf(sysbpm)/dcgain(sysbpm);
    
    [~,fit{bpm}]=compare(m4,detrend(datav,0));
plot(squeeze(fit{bpm})); hold all
bpm
end