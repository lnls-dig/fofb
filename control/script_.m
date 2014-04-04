p = struct('actuator_type', 'current', ...
           'disturbance_amplification', 2, ...
           'phase_margin', 50, ...
           'communication_latency', [80e-6 40e-6 20e-6 10e-6 0], ...
           'ctrl_sampling_rate', [10e3 50e3 100e3], ...
           'corrector_bandwidth', [1e3 2e3 3e3]);

r = fofb_bestctrl(p);

vac_cutoff = 3500;
bpm_group_delay = 4;

best_crossover_frequencies = r.best_crossover_frequencies;

legend_labels = cell(1, size(best_crossover_frequencies, 1)*size(best_crossover_frequencies, 3));
k = 1;
figure;
linestyle = {'-^','--s',':o'};
for i=1:size(best_crossover_frequencies, 3)
    hold on;
    plot(p.communication_latency*1e6, squeeze(best_crossover_frequencies(:,:,i)),linestyle{i});
    [value, unit] = goodunit(p.ctrl_sampling_rate(j), 'Hz', 'tex');
    for j=1:size(best_crossover_frequencies, 3)
        [value2, unit2] = goodunit(p.corrector_bandwidth(j), 'Hz', 'tex');
        legend_labels{k} = sprintf('FOFB data rate = %g %s - Corrector BW = %g %s', value, unit, value2, unit2);
        k = k + 1;
    end
    [value3, unit3] = goodunit(vac_cutoff, 'Hz', 'tex');
    title({sprintf('Vacuum chamber cutoff = %g %s', value3, unit3), sprintf('BPM group delay = %d \\times FOFB samples', bpm_group_delay)});
    xlabel('Latency (\mus)', 'FontSize', 14);
    ylabel('FOFB bandwidth (Hz)', 'FontSize', 14);
    grid on;
end
legend(legend_labels, 'Location', 'Northeast', 'FontSize', 10);
ax = axis;
axis([[min(p.communication_latency) max(p.communication_latency)]*1e6 0 ax(4)]);
set(gca, 'FontSize', 14);