function [sys,fit] = modelo(arquivo, parametros_arx, plano, aquisicao)

%Carrega vetores dos arquivos .mat
A = load(arquivo,'mat_d');
currdata = A.mat_d.currdata;
orbx = double(A.mat_d.orbx);
orby = double(A.mat_d.orby);
sampling_frequency = A.mat_d.sampling_frequency;
prbs_n = A.mat_d.prbs_lfsr_len(1);
prbs_duration = A.mat_d.prbs_step_duration(1);

%Remoção de transiente
currdata = currdata(19051:end, :);
orbx = orbx(19051:end, :);
orby = orby(19051:end, :);

%Separação de conjunto de estimação e validação
h1=length(currdata);
h2=int32(length(currdata)/2);
currdata = currdata(1:h2);
orbx = orbx(1:h2);
orby = orbx(1:h2);

%Vetores para plot de validação
%{
currdata_val = currdata(h2:h1);
orbx_val = orbx(h2:h1);
orby_val = orby(h2:h1);
%}

%Redução do tempo de aquisição (opcional)
T_aq=single(h2)/sampling_frequency;
if exist('aquisicao', 'var')
    k=floor(aquisicao*h2);
    fprintf('Tempo considerado: %f s\n', T_aq*single(k)/single(h2))
    currdata = currdata(1:k);
    orbx = orbx(1:k);
    orby = orby(1:k);
end

%Plot dos sinais
%{
figure;
plot(linspace(0,length(currdata)*sampling_frequency,length(currdata)),currdata)
title('currdata')
figure;
hold on;
plot(linspace(0,length(orbx)*sampling_frequency,length(orbx)),orbx)
plot(linspace(0,length(orby)*sampling_frequency,length(orby)),orby)
title('orbx e orby')
%}

%Filtros (de repositório fofb/sysid)
fir_remove_switching = dfilt.dffir(ones(1,4)/4);
iir_lowpass = design(fdesign.lowpass(0.23, 0.25, 1, 80, 1), 'butter', 'MatchExactly', 'stopband');
currdata = filter(fir_remove_switching, currdata);
orbx = filter(fir_remove_switching, orbx);
orby = filter(fir_remove_switching, orby);
currdata = filter(iir_lowpass, currdata);
orbx = filter(iir_lowpass, orbx);
orby = filter(iir_lowpass, orby);

%Vetores para plot de validação
%{
currdata_val = filter(fir_remove_switching, currdata_val);
orbx_val = filter(fir_remove_switching, orbx_val);
orby_val = filter(fir_remove_switching, orby_val);
currdata_val = filter(iir_lowpass, currdata_val);
orbx_val = filter(iir_lowpass, orbx_val);
orby_val = filter(iir_lowpass, orby_val);
%}

%Retira valor médio
currdata=currdata-mean(currdata);
orbx=orbx-mean(orbx);
orby=orby-mean(orby);

%Vetores para plot de validação
%{
currdata_val=currdata_val-mean(currdata_val);
orbx_val=orbx_val-mean(orbx_val);
orby_val=orby_val-mean(orby_val);
%}

%Subdivide vetor em vetores menores e calcula o valor médio
n=(2^double(prbs_n)-1)*double(prbs_duration); %block size -> siglen (prbs), N*duração do estado

%Ajuste de tamanho
currdata = currdata(1:floor(size(currdata,1)/n)*n, :);
orbx = orbx(1:floor(size(orbx,1)/n)*n, :);
orby = orby(1:floor(size(orby,1)/n)*n, :);

%Vetores para plot de validação
%{
currdata_val = currdata_val(1:floor(size(currdata_val,1)/n)*n, :);
orbx_val = orbx_val(1:floor(size(orbx_val,1)/n)*n, :);
orby_val = orby_val(1:floor(size(orby_val,1)/n)*n, :);
%}

%Obtenção de vetores subdivididos
sublengths = ones(1, floor(numel(currdata)/n))*n;
currdata_div = mat2cell(currdata', 1, sublengths);
aux=zeros(length(sublengths),sublengths(1));
for i=1:numel(sublengths)
    aux(i,:)=currdata_div{1,i};
    currdata_avg=mean(aux,1)';
end
sublengths = ones(1, floor(numel(orbx)/n))*n;
orbx_div = mat2cell(orbx', 1, sublengths);
for i=1:numel(sublengths)
    aux(i,:)=orbx_div{1,i};
    orbx_avg=mean(aux,1)';
end
sublengths = ones(1, floor(numel(orby)/n))*n;
orby_div = mat2cell(orby', 1, sublengths);
for i=1:numel(sublengths)
    aux(i,:)=orby_div{1,i};
    orby_avg=mean(aux,1)';
end

%Vetores para plot de validação
%{
sublengths = ones(1, floor(numel(orbx_val)/n))*n;
orbx_val_div = mat2cell(orbx_val', 1, sublengths);
for i=1:numel(sublengths)
    aux(i,:)=orbx_val_div{1,i};
    orbx_val_avg=mean(aux,1)';
end
sublengths = ones(1, floor(numel(orby_val)/n))*n;
orby_val_div = mat2cell(orby_val', 1, sublengths);
for i=1:numel(sublengths)
    aux(i,:)=orby_val_div{1,i};
    orby_val_avg=mean(aux,1)';
end
sublengths = ones(1, floor(numel(currdata_val)/n))*n;
currdata_val_div = mat2cell(currdata_val', 1, sublengths);
for i=1:numel(sublengths)
    aux(i,:)=currdata_val_div{1,i};
    currdata_val_avg=mean(aux,1)';
end
%}

%Implementação com sigcond (de repositório fofb/sysid)
%{
b = (2^double(prbs_n)-1)*double(prbs_duration); %n blocks -> siglen (prbs), N*duração do estado
currdata_avg = sigcond(currdata, b, iir_lowpass, fir_remove_switching, [0,0]);
orbx_avg = sigcond(orbx, b, iir_lowpass, fir_remove_switching, [0,0]);
%}

%Plot dos sinais avg
%{
figure;
plot(linspace(0,length(currdata_avg)*sampling_frequency*b,length(currdata_avg)),currdata_avg)
title('currdata avg')
figure;
hold on;
plot(linspace(0,length(orbx_avg)*sampling_frequency*b,length(orbx_avg)),orbx_avg)
plot(linspace(0,length(orby_avg)*sampling_frequency*b,length(orby_avg)),orby_avg)
title('orbx e orby avg')
%}

%Plots de visualização
%{
%FFT dos sinais
Fs=sampling_frequency;
L=length(currdata);
Y1=fft(currdata);
figure;
hold on;
plot(Fs/L*(0:L-1),abs(Y1))
L=length(currdata_avg);
Y2=fft(currdata_avg);
plot(Fs/L*(0:L-1),abs(Y2))
title('FFT currdata')

Fs=sampling_frequency;
L=length(orbx);
Y1=fft(orbx);
figure;
hold on;
plot(Fs/L*(0:L-1),abs(Y1))
L=length(orbx_avg);
Y2=fft(orbx_avg);
plot(Fs/L*(0:L-1),abs(Y2))
title('FFT orbx')

Fs=sampling_frequency;
L=length(orby);
Y1=fft(orby);
figure;
hold on;
plot(Fs/L*(0:L-1),abs(Y1))
L=length(orby_avg);
Y2=fft(orby_avg);
plot(Fs/L*(0:L-1),abs(Y2))
title('FFT orby')
%}%{
%Periodograma currdata
Fs=sampling_frequency;
L=length(currdata);
[Pxx,F] = periodogram(currdata,hamming(length(currdata)),length(currdata),Fs);
figure;
hold on;
plot(F,10*log10(Pxx))

Fs=sampling_frequency;
L=length(currdata_avg);
[Pxx,F] = periodogram(currdata_avg,hamming(length(currdata_avg)),length(currdata_avg),Fs);
%figure;
plot(F,10*log10(Pxx))
title('Periodograma currdata')

%Periodograma orbx
Fs=sampling_frequency;
L=length(orbx);
[Pxx,F] = periodogram(orbx,hamming(length(orbx)),length(orbx),Fs);
figure;
hold on;
plot(F,10*log10(Pxx))

Fs=sampling_frequency;
L=length(orbx_avg);
[Pxx,F] = periodogram(orbx_avg,hamming(length(orbx_avg)),length(orbx_avg),Fs);
%figure;
plot(F,10*log10(Pxx))
title('Periodograma orbx')

%Periodograma orby
Fs=sampling_frequency;
L=length(orby_H);
[Pxx,F] = periodogram(orby,hamming(length(orby)),length(orby),Fs);
figure;
hold on;
plot(F,10*log10(Pxx))

Fs=sampling_frequency;
L=length(orby_H_avg);
[Pxx,F] = periodogram(orby_avg,hamming(length(orby_avg)),length(orby_avg),Fs);
%figure;
plot(F,10*log10(Pxx))
title('Periodograma orby')
%}

%Aplicação no modelo ARX
if plano=='X'
    sys = arx(currdata_avg,orbx_avg,parametros_arx, arxOptions('Focus', 'Simulation', 'EnforceStability', true));
elseif plano=='Y'
    sys = arx(currdata_avg,orby_avg,parametros_arx, arxOptions('Focus', 'Simulation', 'EnforceStability', true));
else 
    print('Erro no parâmetro plano')
end
fit=sys.Report.Fit.FitPercent;

%Plots básicos
%{
[y,tOut] = step(sys);
figure;
plot(tOut,y)
figure;
bode(sys)
figure;
compare([z1(5:end,1) z1(5:end,2)],sys,Inf)
%}

%Plot do sistema e conjunto de validação
%{
Ts=T_aq*aquisicao*2;
t=linspace(0,Ts/2,length(double(currdata_avg)))';
y = lsim(sys, double(currdata_avg), []);
figure;
plot(t,y,t,orbx_avg)
title('Estimação')
t=linspace(Ts/2,Ts,length(double(currdata_val_avg)))';
y = lsim(sys, double(currdata_val_avg), []);
figure;
plot(t,y,t,orbx_val_avg)
title('Extrapolação')
%}