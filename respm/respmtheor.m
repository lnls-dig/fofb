function M = respmtheor(beta, mu, eta, alpha, L, tune_t, orbit_indexes, corr_indexes, type, tau_t)

if nargin < 9 || isempty(type)
    type = 'static';
end

if nargin < 10 || isempty(tau_t)
    % Transverse plane damping time [turns]
    tau_t = Inf;
end

% Relative phase advances among each orbit and corrector (or dipolar field) points
ph_orbit = repmat(mu(orbit_indexes), 1, length(corr_indexes));
ph_corr = repmat(mu(corr_indexes)', length(orbit_indexes), 1);
ph_orbit_minus_corr = ph_orbit-ph_corr;

% Dispersive orbit term
eta_orbit = repmat(eta(orbit_indexes), 1, length(corr_indexes));
eta_corr = repmat(eta(corr_indexes)', length(orbit_indexes), 1);
dispersive_orbit = eta_orbit.*eta_corr/alpha/L;

% Build response matrix

if strcmpi(type, 'static')
    % Betatron phase advance term
    M = cos(pi*tune_t - abs(ph_orbit_minus_corr));
    
    % Tune term
    M = M/2/sin(pi*tune_t);
elseif strcmpi(type, 'dynamic')
    put_delay = ph_orbit_minus_corr < 0;
    ph_orbit_minus_corr(ph_orbit_minus_corr < 0) = ph_orbit_minus_corr(ph_orbit_minus_corr < 0) + 2*pi*tune_t;
    
    q = exp(-1/tau_t);
    ct = cos(2*pi*tune_t);
    st = sin(2*pi*tune_t);
    
    den = [1 -2*q*ct q^2];
    num_sin = [0 q*st 0];
    num_cos = [1 -q*ct 0];

    C = cos(ph_orbit_minus_corr);
    S = sin(ph_orbit_minus_corr);
    Cdelay = zeros(size(C));
    Sdelay = zeros(size(S));
    Cdelay(put_delay) = C(put_delay);
    Sdelay(put_delay) = S(put_delay);
    C(put_delay) = 0;
    S(put_delay) = 0;

    num_ = zeros(size(ph_orbit_minus_corr, 1), 3*size(ph_orbit_minus_corr, 2));
    num_(:, 1:3:end) = S;
    num_(:, 2:3:end) = C*num_sin(2) + S*num_cos(2) + Sdelay;
    num_(:, 3:3:end) = Cdelay*num_sin(2) + Sdelay*num_cos(2);
    num = mat2cell(num_, ones(size(ph_orbit_minus_corr, 1),1), 3*ones(size(ph_orbit_minus_corr, 2),1));
    M = tf(num, den, 1);
else
    error('Non-valid ''type'' argument.');
end

% Betatron "modulation" term
sqrt_beta_orbit = sqrt(beta(orbit_indexes));
sqrt_beta_cm = sqrt(beta(corr_indexes));

for i=1:length(orbit_indexes)
    M(i,:) = M(i,:).*sqrt_beta_cm';
end
for j=1:length(corr_indexes)
    M(:,j) = M(:,j).*sqrt_beta_orbit;
end

% Dispersive orbit term and additional treatments
if strcmpi(type, 'static')
    M = M + dispersive_orbit;
elseif strcmpi(type, 'dynamic')
    M = ss(M);
    M = balred(M, 2);
    
    % TODO: include dispersive orbit dynamic response
    M = M + dispersive_orbit;
end