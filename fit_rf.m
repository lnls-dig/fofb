% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
%Input data from experiment
y_data = h5read('rf-phase-fofb-off.h5','/data/orbx');
sampling_frequency = h5read('rf-phase-fofb-off.h5','/data/sampling_frequency');
Ts=1/sampling_frequency;
fir_moving_average = dfilt.dffir(ones(1,4)/4);
%Defines model for fitting
X=@(x,xdata) x(1).*exp(-x(2).*xdata).*sin(x(3).*xdata)+x(4).*exp(-x(2).*xdata).*cos(x(3).*xdata);
x0=[100 0 2*pi*2000 100];
alpha_s = 0;
omega = 0;
%Attempts 15 fits
%BPM 19, step inputs chosen visually, interval manually chosen
for i=1:15
    %Positive
    if i==1
        y = double(y_data(19,13081:13197)); 
    elseif i==2
        y = double(y_data(19,13000+2164:13000+2361));
    elseif i==3
        y = double(y_data(19,13000+4710:13000+5241));
    elseif i==4
        y = double(y_data(19,13000+9338:13000+9573));
    elseif i==5
        y = double(y_data(19,13000+18602:13000+18821));
    elseif i==6
        y = double(y_data(19,13000+106562:13000+106745));
    elseif i==7
        y = double(y_data(19,13000+111190:13000+111329));
    elseif i==8
        y = double(y_data(19,13000+115822:13000+115989));
    elseif i==9
        y = double(y_data(19,13000+120450:13000+120601));
    elseif i==10
        y = double(y_data(19,13000+152858:13000+153001));
    %Negative
    elseif i==11
        y = double(y_data(19,13000+154941:13000+155069));
    elseif i==12
        y = double(y_data(19,13000+159572:13000+159721));
    elseif i==13
        y = double(y_data(19,13000+164201:13000+164389));
    elseif i==14
        y = double(y_data(19,13000+168828:13000+168891));
    elseif i==15
        y = double(y_data(19,13000+173460:13000+173629));
    end
%y = filter(fir_moving_average, y); %optional filter
y = detrend(y,0);
%y=y-mean(y);
t = linspace(0,length(y)*Ts,length(y));
%Performs fit
[x,resnorm,~,exitflag,output] = lsqcurvefit(X,x0,t,y);
figure;
plot(t,X(x,t),t,y);
grid on
fprintf("Step %d\n", i);
fprintf("Omega: %f Hz\n", x(3)/(2*pi));
fprintf("alpha_s: %f s^-1\n", x(2));
alpha_s=alpha_s+x(2);
omega=omega+x(3);
end
alpha_s=alpha_s/15;
omega=omega/15;
omega_hz = omega/2/pi;
fitted_parameters.omega = omega;
fitted_parameters.omega_hz = omega_hz;
fitted_parameters.alpha_s = alpha_s;
save('rf_fitted_parameters','fitted_parameters');