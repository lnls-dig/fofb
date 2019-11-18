if ~exist('THERING','var')
    global THERING;
    sirius;
    setoperationalmode(1);
end

%%
coupling = 1/100;

unitconv = 1e-6;

fprintf('   Calculating beam size\n'); tic;
twiss = calctwiss;
eqp = atsummary; 
betax = twiss.betax;
betay = twiss.betay;
etax = twiss.etax;
sigmax = sqrt(eqp.naturalEmittance*betax + (eqp.naturalEnergySpread*etax).^2)/unitconv;
sigmay = sqrt(coupling*eqp.naturalEmittance*betay)/unitconv;
sigmax_ = sqrt(eqp.naturalEmittance./betax + (eqp.naturalEnergySpread*etax).^2)/unitconv;
sigmay_ = sqrt(coupling*eqp.naturalEmittance./betay)/unitconv;
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);


mia  = findcells(THERING, 'FamName', 'mia');  % straight section A
mib  = findcells(THERING, 'FamName', 'mib');  % straight section B
mip  = findcells(THERING, 'FamName', 'mip');  % straight section P
mc  = findcells(THERING, 'FamName', 'mc');    % bending magnet center
bpm = findcells(THERING, 'FamName', 'BPM');  % BPMs
    
light_sources =  sort([mc mia mib mip]);

family_data = sirius_si_family_data(THERING);
hcm_slow = sort(family_data.CH.ATIndex);
vcm_slow = sort(family_data.CV.ATIndex);
hcm_fast = sort([family_data.FC1.ATIndex; family_data.FC2.ATIndex]);
vcm_fast = hcm_fast;



th_inf = mia(2);
th_sup = mib(3);
all_points = th_inf:th_sup;
light_sources = light_sources(light_sources >= th_inf & light_sources <= th_sup);
all_points = all_points(all_points >= th_inf & all_points <= th_sup);
hcm_slow = hcm_slow(hcm_slow >= th_inf & hcm_slow <= th_sup);
vcm_slow = vcm_slow(vcm_slow >= th_inf & vcm_slow <= th_sup);
hcm_fast = hcm_fast(hcm_fast >= th_inf & hcm_fast <= th_sup);
vcm_fast = vcm_fast(vcm_fast >= th_inf & vcm_fast <= th_sup);


pos0 = twiss.pos(th_inf);
%%
figure;
semilogy(twiss.pos(all_points)-pos0, sigmax(all_points), 'b', 'LineWidth', 1.2);
hold on;
semilogy(twiss.pos(all_points)-pos0, sigmay(all_points), 'r', 'LineWidth', 1.2);
semilogy(twiss.pos(light_sources)-pos0, sigmax(light_sources),'xb','MarkerSize',20);
semilogy(twiss.pos(light_sources)-pos0, sigmay(light_sources),'xr','MarkerSize',20);
semilogy(twiss.pos(bpm)-pos0,  sigmax(bpm),'o','MarkerEdgeColor','m','MarkerFaceColor','m');
semilogy(twiss.pos(hcm_slow)-pos0, sigmax(hcm_slow),'s','MarkerEdgeColor','b','MarkerFaceColor',[0.6 0.6 1],'MarkerSize',8);
semilogy(twiss.pos(vcm_slow)-pos0, sigmay(vcm_slow),'s','MarkerEdgeColor','r','MarkerFaceColor',[1 0.6 0.6],'MarkerSize',8);
semilogy(twiss.pos(hcm_fast)-pos0, sigmax(hcm_fast),'p','MarkerEdgeColor',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',10);
semilogy(twiss.pos(bpm)-pos0,  sigmay(bpm),'o','MarkerEdgeColor','m','MarkerFaceColor','m');
semilogy(twiss.pos(hcm_fast)-pos0, sigmay(hcm_fast),'p','MarkerEdgeColor',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',10);
legend('Beam Size (H plane)', 'Beam Size (V plane)', 'Light Sources (H)', 'Light Sources (V)', 'BPMs', 'Slow Correctors (H)', 'Slow Correctors (V)', 'Fast Correctors');

axis([twiss.pos(all_points(1))-1-pos0 twiss.pos(all_points(end))+1-pos0 [1e-6 1e-4]/unitconv]);
xlabel('s [m]')
ylabel('Beam size [\mum]')
grid on
%
%%
figure;
semilogy(twiss.pos(all_points)-pos0, sigmax_(all_points), 'b', 'LineWidth', 1.2);
hold on;
semilogy(twiss.pos(all_points)-pos0, sigmay_(all_points), 'r', 'LineWidth', 1.2);
semilogy(twiss.pos(light_sources)-pos0, sigmax_(light_sources),'xb','MarkerSize',20);
semilogy(twiss.pos(light_sources)-pos0, sigmay_(light_sources),'xr','MarkerSize',20);
semilogy(twiss.pos(bpm)-pos0,  sigmax_(bpm),'o','MarkerEdgeColor','m','MarkerFaceColor','m');
semilogy(twiss.pos(bpm)-pos0,  sigmay_(bpm),'o','MarkerEdgeColor','m','MarkerFaceColor','m');


% figure;
% semilogy(twiss.pos(all_points), [sigmax_(all_points) sigmay_(all_points)]);
% hold on; semilogy(twiss.pos(light_sources),  [sigmax_(light_sources) sigmay_(light_sources)],'.');
% hold on; semilogy(twiss.pos(bpm),  [sigmax_(bpm) sigmay_(bpm)],'ok');
axis([twiss.pos(all_points(1))-1-pos0 twiss.pos(all_points(end))+1-pos0 [1e-7 1e-4]/unitconv]);
xlabel('s [m]')
ylabel('Beam divergence [\murad]')
grid on
legend('Beam Divergence (H plane)', 'Beam Divergence (V plane)', 'Light Sources (H)', 'Light Sources (V)', 'BPMs');