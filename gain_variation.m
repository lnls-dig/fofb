% Copyright (C) 2024 CNPEM (cnpem.br)
% Author: Lucas Pelike <lucas.pelike@lnls.br>

clc;clear;
%Gain  variation
%Simulation parameters
min_gain = 0.001;
gain_step = 0.025;
max_gain = 0.5;
%FOFB parameters
load respmat.mat
M = mat_d';
sys=cell(size(M,2), 1);
load sysid_res-24-08-13.mat
Ts=1/48193;
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
fprintf("Elapsed time: %f s\n", toc);
%Considering correctors discrepancies
%stable = isstable(P);
it_h=1;
stable_complete = zeros(int32(max_gain/gain_step));
for i=min_gain:gain_step:max_gain
    it_v=1;
    fprintf("It_h = %d \n", it_h);
    for j=min_gain:gain_step:max_gain
        fprintf("It_v = %d \n", it_v);
        K = new_gain(i,j,size(M,2)/2);
        fprintf("K_h = %f , K_v = %f\n", K(1)/48193, K(81)/48193);
        [P, ~, ~] = ofbmdl(M, Mc, K, sys);
        stable_complete(it_h,it_v) = isstable(P);
        it_v=it_v+1;
    end
    it_h=it_h+1;
end
%Equal correctors
sys = sys(10); %Use corrector 10 as all correctors
it_h=1;
stable_simple = zeros(int32(max_gain/gain_step));
for i=min_gain:gain_step:max_gain
    it_v=1;
    fprintf("It_h = %d \n", it_h);
    for j=min_gain:gain_step:max_gain
        fprintf("It_v = %d \n", it_v);
        K = new_gain(i,j,size(M,2)/2);
        fprintf("K_h = %f , K_v = %f\n", K(1)/48193, K(81)/48193);
        [P, ~, ~] = ofbmdl(M, Mc, K, sys);
        stable_simple(it_h,it_v) = isstable(P);
        it_v=it_v+1;
    end
    it_h=it_h+1;
end

%Plots results
%Dots plots
figure;
title('Stability plot for different horizontal and vertical gains')
hold on;
it_h=1;
for i=min_gain:gain_step:max_gain
    it_v=1;
    for j=min_gain:gain_step:max_gain
        if stable_complete(it_h,it_v)==1
            plot(i,j,'Color','red','Marker','.')
        else
            plot(i,j,'Color','black','Marker','.')
        end
        if stable_simple(it_h,it_v)==1
            plot(i,j,'Color','blue','Marker','d')
        else
            plot(i,j,'Color','black','Marker','d')
        end
        it_v=it_v+1;
    end
    it_h=it_h+1;
end
grid on;
xlabel('Vertical plane gain')
ylabel('Horizontal plane gain')

%Area plots
figure;
title('Stability plot for different horizontal and vertical gains')
hold on;
idx_complete_v = find(stable_complete(1,:) == 0,1,'first');
idx_complete_h = find(stable_complete(:,1) == 0,1,'first');
idx_simple_v = find(stable_simple(1,:) == 0,1,'first');
idx_simple_h = find(stable_simple(:,1) == 0,1,'first');
it_h=1;
it_v=1;
for i=min_gain:gain_step:max_gain
    it_v=1;
    for j=min_gain:gain_step:max_gain
        if it_v==idx_complete_v
            lim_complete_v=j;
        end
        if it_h==idx_complete_h
            lim_complete_h=i;
        end
        if it_v==idx_simple_v
            lim_simple_v=j;
        end
        if it_h==idx_simple_h
            lim_simple_h=i;
        end
        it_v=it_v+1;
    end
    it_h=it_h+1;
end
area([0 lim_simple_v],[lim_simple_h lim_simple_h],'FaceColor','red','FaceAlpha',0.6);
area([0 lim_complete_v],[lim_complete_h lim_complete_h],'FaceColor','blue');
grid on;
xlabel('Vertical plane gain')
ylabel('Horizontal plane gain')
legend({'Matched actuators','Discrepant actuators'})

function [K] = new_gain(h_gain, v_gain, plane_size)
    K = [ones(plane_size,1)*h_gain*48193; ones(plane_size,1)*v_gain*48193];
end