function M = respmtheor(beta, mu, eta, alpha, L, tune, orbit_indexes, corr_indexes)

% Relative phase advances among each orbit and corrector (or dipolar field) points
ph_orbit = repmat(mu(orbit_indexes), 1, length(corr_indexes));
ph_corr = repmat(mu(corr_indexes)', length(orbit_indexes), 1);
ph_orbit_minus_corr = ph_orbit-ph_corr;

% Dispersive orbit term
eta_orbit = repmat(eta(orbit_indexes), 1, length(corr_indexes));
eta_corr = repmat(eta(corr_indexes)', length(orbit_indexes), 1);
dispersive_orbit = eta_orbit.*eta_corr/alpha/L;

% Build response matrix

% Betatron phase advance term
M = cos(pi*tune - abs(ph_orbit_minus_corr));

% Tune term
M = M/2/sin(pi*tune);

% Betatron "modulation" term
sqrt_beta_orbit = sqrt(beta(orbit_indexes));
sqrt_beta_cm = sqrt(beta(corr_indexes));

for i=1:length(orbit_indexes)
    M(i,:) = M(i,:).*sqrt_beta_cm';
end
for j=1:length(corr_indexes)
    M(:,j) = M(:,j).*sqrt_beta_orbit;
end

% Dispersive orbit term
M = M + dispersive_orbit;