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
sigmax = sqrt(eqp.naturalEmittance*betax + (eqp.naturalEnergySpread*twiss.etax).^2)/unitconv;
sigmay = sqrt(coupling*eqp.naturalEmittance*betay)/unitconv;
sigmax_ = sqrt(eqp.naturalEmittance./betax + (eqp.naturalEnergySpread*twiss.etax).^2)/unitconv;
sigmay_ = sqrt(coupling*eqp.naturalEmittance./betay)/unitconv;
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);


mia  = findcells(THERING, 'FamName', 'mia');  % straight section A
mib  = findcells(THERING, 'FamName', 'mib');  % straight section B
mip  = findcells(THERING, 'FamName', 'mip');  % straight section P
mc  = findcells(THERING, 'FamName', 'mc');    % bending magnet center
bpm  = findcells(THERING, 'FamName', 'BPM');    % BPMs

light_source_points =  sort([mc mia mib mip]);
all_points = sort([bpm light_source_points]);

th_inf = mia(1);
th_sup = mia(2);
light_source_points = light_source_points(light_source_points >= th_inf & light_source_points <= th_sup);
all_points = all_points(all_points >= th_inf & all_points <= th_sup);

figure;
semilogy(twiss.pos(all_points), [sigmax(all_points) sigmay(all_points)], '.-');
hold on; semilogy(twiss.pos(light_source_points),  [sigmax(light_source_points) sigmay(light_source_points)],'o');
axis([twiss.pos(all_points(1))-1 twiss.pos(all_points(end))+1 [1e-6 1e-4]/unitconv]);
xlabel('s [m]')
ylabel('Beam size [\mum]')
grid on
legend('H plane', 'V plane', 'H plane (only IDs)', 'V plane (only IDs)');

figure;
semilogy(twiss.pos(all_points), [sigmax_(all_points) sigmay_(all_points)], '.-');
hold on; semilogy(twiss.pos(light_source_points),  [sigmax_(light_source_points) sigmay_(light_source_points)],'o');
axis([twiss.pos(all_points(1))-1 twiss.pos(all_points(end))+1 [1e-7 1e-4]/unitconv]);
xlabel('s [m]')
ylabel('Beam divergence [\murad]')
grid on
legend('H plane', 'V plane', 'H plane (only light sources)', 'V plane (only light sources)');