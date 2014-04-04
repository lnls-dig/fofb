if ~exist('THERING','var')
global THERING;
sirius;
setoperationalmode(1);
end


%%
coupling = 1/100;

fprintf('   Calculating beam size\n'); tic;
twiss = calctwiss;
eqp = atsummary; 
betax = twiss.betax;
betay = twiss.betay;
sigmax = sqrt(eqp.naturalEmittance*betax + (eqp.naturalEnergySpread*twiss.etax).^2);
sigmay = sqrt(coupling*eqp.naturalEmittance*betay);
fprintf('   Done. Elapsed time: %0.3f s\n\n', toc);


mia  = findcells(THERING, 'FamName', 'mia');  % straight section A
mib  = findcells(THERING, 'FamName', 'mib');  % straight section B
mc  = findcells(THERING, 'FamName', 'mc');  % bending magnet center
id_points = sort([mia mib]);
light_source_points =  sort([mc id_points]);

figure;
plot(light_source_points, [sigmax(light_source_points) sigmay(light_source_points)]);
hold on; plot(id_points,  [sigmax(id_points) sigmay(id_points)],'o');