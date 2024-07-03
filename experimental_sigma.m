%FOFB model with PRBS inputs

clc;clear;
%Simulation parameters
step_duration = 2;
n_prbs=8;
prbs_periods = 5; %Standard is 24
prbs_amplitude = 4000; %Standard is 4000;
standard_deviation = 300;
type='Open Loop'; %Sensibility or Open Loop
processing_mode='Averaging'; %Averaging or PRBS_binning
%FOFB parameters
load respmat.mat
load sysid_res.mat
M = mat_d';
Ts=1/48193;
Ts_prbs=step_duration*Ts;
K = 0.052*48193;
fofb_type.bpm_sector = 'm1m2c2c3';
fofb_type.corr_sector = 'm1m2c2c3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = [1 80 81 160];
[bpm_idx, corr_idx] = fofb_idx(fofb_type);
%Normalizing models and removing 'NaN' values
fprintf("Normalizing sysid models...\n");
tic
for i=1:size(M,2)
    if (i==1)||(i==80)||(i==81)||(i==160) %excluded correctors: 1, 80, 81 and 160
        sys{i} = tf(0,1,'Ts',Ts);
    else
        sys{i} = idtf(sys{i}/dcgain(idtf(sys{i})),'Ts',Ts);
    end
end
fprintf("Elapsed time: %f s\n", toc);
fprintf("Building model...\n");
tic
M = M(bpm_idx, corr_idx);
Mc = pinv(M);
sys = sys(corr_idx);
[P, G, C] = ofbmdl(M, Mc, K, sys);
fprintf("Elapsed time: %f s\n", toc);
%PRBS input
fprintf("Creating PRBS input...\n");
N=2^n_prbs-1;
fprintf("Acquisition time (mode): %f s\n", N*prbs_periods*Ts);
input = frest.PRBS('Order',n_prbs,'NumPeriods',prbs_periods,'Amplitude',prbs_amplitude,'Ts',Ts,'UseWindow','off');
inputts = generateTimeseries(input);
prbs_input = repelem(inputts.Data,step_duration);
prbs_input = prbs_input-mean(prbs_input);
prbs_time=0:Ts:length(prbs_input)*Ts;
prbs_time(end)=[];
plot_signal(prbs_input,Ts);
fprintf("Creating input vectors...\n");
%Defining U, V, u and y
if strcmp(type,'Sensibility')
    %SYS = P('yd','d');
    SYS=P('yd',{'d','n'});
    %[U,Sigma,V] = svd(dcgain(SYS));
    [U,Sigma,V] = svd(M);
    n_excited_modes =  min(size(U,1),size(V,1));
    m = size(U,1);
    n = m;
    u=zeros(length(prbs_time),n,n_excited_modes);
    y=zeros(length(prbs_time),m,n_excited_modes);
    tic
    for i=1:n_excited_modes
        u(:,:,i) = U(:,i)'.*prbs_input;
    end
    fprintf("Elapsed time: %f s\n", toc);
elseif strcmp(type,'Open Loop')
    SYS = G;
    %[U,Sigma,V] = svd(dcgain(SYS));
    [U,Sigma,V] = svd(M);
    n_excited_modes =  min(size(U,1),size(V,1));
    m = size(U,1);
    n = size(V,1);
    u=zeros(length(prbs_time),n,n_excited_modes);
    y=zeros(length(prbs_time),m,n_excited_modes);
    tic
    for i=1:n_excited_modes
        u(:,:,i) = V(:,i)'.*prbs_input;
    end
    fprintf("Elapsed time: %f s\n", toc);
end
%Plot sample input signal
plot_signal(u(:,1,1),Ts);
%Simulating system
tic
%Define frequency indexes for PRBS binning
if strcmp(processing_mode,'PRBS_binning')
    [~,F_index] = findpeaks(abs(adj_fft(u(:,1,1),Ts)),'MinPeakHeight',0.1e-10); %Finds PRBS peaks
elseif strcmp(processing_mode,'Averaging')
    u_avg = reshape(u(:,1,1),N*step_duration,prbs_periods);
    u_avg = mean(u_avg(:,2:end),2);
    [~,F] = adj_fft(u_avg,Ts);
    F_index = 1:1:length(F);
end
fprintf("Simulating system and creating FFT of input and output vectors (Y_u(f) and Y_y(f))... \n");
Yu = ones(length(F_index),n,n_excited_modes);
Yy = ones(length(F_index),m,n_excited_modes);
Yu_aux = zeros(n,n_excited_modes);
Yy_aux = zeros(m,n_excited_modes);
snr_db = ones(m,n_excited_modes);
%Iterates for each mode tested
for k=1:n_excited_modes
    fprintf("Mode:%d\n",k);
    if strcmp(type,'Sensibility')
        noise = standard_deviation.*randn(length(y(:,1,1)),length(y(1,:,1)));
        y(:,:,k) = lsim(SYS,[u(:,:,k) noise],prbs_time);
    elseif strcmp(type,'Open Loop')
        y(:,:,k) = lsim(SYS,u(:,:,k),prbs_time);
        noise = standard_deviation.*randn(length(y(:,1,1)),1);
        y(:,i,k) = y(:,i,k) + noise;
    end
    %Iterates in each input/output
    for i=1:m
        %Building Yu arrays
        if(i<=n)
            if strcmp(processing_mode,'PRBS_binning')
                [Yu_aux,~] = adj_fft(u(:,i,k),Ts);
            elseif strcmp(processing_mode,'Averaging')
                u_avg = reshape(u(:,i,k),N*step_duration,prbs_periods);
                u_avg = mean(u_avg(:,2:end),2);
                [Yu_aux,~] = adj_fft(u_avg,Ts);
            end
            Yu(:,i,k) = Yu_aux(F_index);
        end
        %Building Yy arrays
        if strcmp(processing_mode,'PRBS_binning')
            [Yy_aux,F] = adj_fft(y(:,i,k),Ts);   
            Yy(:,i,k) = Yy_aux(F_index);
        elseif strcmp(processing_mode,'Averaging')
            y_avg = reshape(y(:,i,k),N*step_duration,prbs_periods);
            y_avg = mean(y_avg(:,2:end),2);
            [Yy_aux,F] = adj_fft(y_avg,Ts);   
            Yy(:,i,k) = Yy_aux(F_index);
        end
        %Calculating SNR vector
        if strcmp(type,'Open Loop')
            snr_db(i,k) = 20*log10(rms(y(:,i,k))/rms(noise));
        elseif strcmp(type,'Sensibility')
            snr_db(i,k) = 20*log10(rms(y_avg)/rms(noise(:,i)));
        end
        %Print value for first output
        if i==1
            fprintf("Output 1 signal standard deviation: %f\n", std(detrend(y(:,1,k),0)));
            if strcmp(processing_mode,'Averaging')
                fprintf("Output 1 signal standard deviation after averaging: %f\n", std(detrend(y_avg,0)));
            end
        end
    end 
end
%Plot SNR for some modes
figure;
hold on;
s = 1:1:m;
for k=1:20:n_excited_modes
    plot(s,snr_db(:,k))
end
title("SNR for subset of tested modes, "+ prbs_periods +" PRBS periods, " +"N PRBS = "+n_prbs)
xlabel('Output')
ylabel('SNR (dB)')
%Plot sample output signal
plot_signal(y(:,1,1),Ts);
%Clear unused variables to clear space
clear Yu_aux Yy_aux U V u y prbs_input prbs_time;
%Calculating pseudo-inverses and SVDs
div = ones(size(Yu,2),size(Yu,3));
aux = ones(size(Yy,2),size(Yy,3));
sigma_exp = ones(m,n,length(F_index));
center_f = find(abs(Yu(:,1,1))==max(abs(Yu(:,1,1)))); 
for f=center_f(2):length(F_index)
    div(:,:) = Yu(f,:,:);
    aux(:,:) = Yy(f,:,:);
    [~,s_f,~] = svd(aux*pinv(div));
    sigma_exp(:,:,f)=s_f;
end
fprintf("Elapsed time: %f s\n", toc);
%Plot results
if strcmp(type,'Sensibility')
    SYS = P('yd','d');
elseif strcmp(type,'Open Loop')
    SYS = G;
end
figure;
subplot(1,3,1);
hold on;
sigma_aux = ones(min(m,n),length(F_index));
for k=1:min(m,n)
    sigma_aux(k,:) = abs(sigma_exp(k,k,:)); 
    plot(F(F_index(center_f(2):end)),sigma_aux(k,center_f(2):end));
end
title('Singular Values - linear')
xlabel('Hz')
ylabel('Singular value')
hold on;
for k=1:min(m,n)
    sigma_aux(k,:) = abs(sigma_exp(k,k,:));
    semilogx(2*pi*F(F_index(center_f(2):end)),20*log10(sigma_aux(k,center_f(2):end)));
end
title('Singular Values - logx20*log')
xlabel('rad/s')
ylabel('Singular value (dB)')
ax = gca;
ax.XScale ='log';
subplot(1,3,2);
hold on;
for k=1:min(m,n)
    sigma_aux(k,:) = abs(sigma_exp(k,k,:));
    plot(F(F_index(center_f(2):end)),log10(sigma_aux(k,center_f(2):end)));
end
title('Singular Values - semilog')
xlabel('Hz')
ylabel('Singular value (dB)')
subplot(1,3,3);
sigma(SYS);
figure;
hold on;
for k=1:min(m,n)
    sigma_aux(k,:) = abs(sigma_exp(k,k,:));
    semilogx(2*pi*F(F_index(center_f(2):end)),20*log10(sigma_aux(k,center_f(2):end)),"green");
end
title('Singular Values - Experimental (green) x Sigma (blue)')
xlabel('rad/s')
ylabel('Singular value (dB)')
ax = gca;
ax.XScale ='log';
hold on;
sigma(SYS);

function [] = plot_signal(x, Ts)
figure;
Fs=1/Ts;
L=length(x);
%Periodogram
subplot(1,3,1);
[Pxx,F] = periodogram(x,hamming(L),L,Fs);
plot(F,10*log10(Pxx))
title('Periodogram')
xlabel('Hz')
ylabel('dB')
%FFT 
subplot(1,3,2);
Y=fft(x);                       
Y=fftshift(Y);                  
Y=abs(Y);                       
Y=Y/L;                         
F=(((-L/2:L/2-1))')*(Fs/L);     
plot(F,Y)
title('FFT')
xlabel('Hz')
%After PRBS binning
subplot(1,3,3);
[~,F_index] = findpeaks(Y,'MinPeakHeight',0.1e-10);
plot(F(F_index),Y(F_index));
title('FFT - separated PRBS bins')
xlabel('Hz')
end

function [Y, F] = adj_fft(x,Ts)
Fs=1/Ts;
L=length(x);
Y=fft(x);                       
Y=fftshift(Y);                  
%Y=abs(Y);                       
%Y=Y/L;                         
F=(((-L/2:L/2-1))')*(Fs/L);
end