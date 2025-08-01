function [F, AF, nf] = eqfilt(A, pz_left, pz_right, znmp_left, tgt_bw, max_bw, n)
% EQFILT Equalization filters.
% 
% Find discrete-time equalization filters to achieve a target bandwidth.
% Make poles and zeros cancellations inside a region of interest and add
% one missing non-minimum phase zero close to -1 if needed.
%
% [F, AF, nf] = eqfilt(A, pz_left, pz_right, znmp_left, tgt_bw, max_bw, n)
%
% INPUTS:
%   A:         Open loop transfer functions corrector by corrector (cell
%              array of dynamical system objects).
%   pz_left:   Left edge of region of interest in the pole-zero map for
%              pole-zero cancellations (default: -1).
%   pz_right:  Right edge of region of interest in the pole-zero map for
%              pole-zero cancellations (default: 1).
%   znmp_left: Left edge of non-minimum phase zeros region of interest in
%              the pole-zero map (default: -20).
%   tgt_bw:    Target banwdidth (dominant pole, -3 dB cutoff) [Hz]
%              (default: inf).
%   max_bw:    Banwdidth of poles which need to be added to make filters
%              realizable after addition of zeros (default: inf).
%   n:         Target filters order. The resulting filter order is reduced
%              with BALRED when necessary.
%
% OUTPUTS:
%   F:    Resulting equalization filters.
%   AF:   Combined response (open loop responses and filters).
%   nf:   Filters orders before order reduction.

if nargin < 2 || isempty(pz_left) || pz_left < 0
    pz_left = -1;
end

if nargin < 3 || isempty(pz_right) || pz_right < 0
    pz_right = 1;
end

if nargin < 4 || isempty(znmp_left) || znmp_left < 0
    znmp_left = -20;
end

if nargin < 5 || isempty(tgt_bw) || tgt_bw < 0
    tgt_bw = inf;
end

if nargin < 6 || isempty(max_bw) || max_bw < 0
    max_bw = inf;
end

if nargin < 7 || isempty(n) || n < 0
    n = 0;
end

% Detrmine maximum number of zeros and poles in transfer functions of A
nz = 0;
np = 0;
Azpk = cell(size(A));
for i=1:length(A)
    Azpk{i} = zpk(minreal(A{i}));
    nz_ = length(Azpk{i}.z{1});
    np_ = length(Azpk{i}.p{1});
    if nz_ > nz
        nz = nz_;
    end
    if np_ > np
        np = np_;
    end
end

% Create matrices used to manipulate poles and zeros
p = nan(length(A),np);          % Poles
z = nan(length(A),nz);          % Zeros
p_rem = nan(length(A),np);      % Poles to be canceled by adding zeros in the filter F
z_rem = nan(length(A),nz);      % Zeros to be canceled by adding zeros in the filter F
nmp_z_add = nan(length(A),nz);  % Non-minimum phase zeros to be added

% Extract all poles and zeros of the actuator responses A to be equalized
for i=1:length(A)
    z_ = Azpk{i}.z{1};
    p_ = Azpk{i}.p{1};
    z(i,1:length(z_)) = z_;
    p(i,1:length(p_)) = p_;
end

% Find all poles and zeros to be canceled inside the region of interest
z_rem(real(z) >= pz_left & real(z) <= pz_right) = z(real(z) >= pz_left & real(z) <= pz_right);
p_rem(real(p) >= pz_left & real(p) <= pz_right) = p(real(p) >= pz_left & real(p) <= pz_right);

% Find dominant non-minimum phase zeros
nmp_z_add(real(z) < -1 & real(z) >= znmp_left) = z(real(z) < -1 & real(z) >= znmp_left);
idx_add_z = find(all(isnan(nmp_z_add),2));

% Compute average value of non-minimum phase zeros
nmp_z_add = nmp_z_add(~isnan(nmp_z_add));
nmp_z_add = mean(nmp_z_add);

nf = zeros(length(A),1);
F = cell(size(A));
AF = cell(size(A));
for i=1:length(A)
    Ts = A{i}.Ts;
    
    % Isolate poles and zeros to be canceled
    z_ = z_rem(i,:);
    z_ = z_(~isnan(z_));
    p_ = p_rem(i,:);
    p_ = p_(~isnan(p_));
    
    % Add non-minimum phase zeros that may be missing
    if ismember(i, idx_add_z)
        p_ = [p_ nmp_z_add];
    end
    
    % Replace poles to the target bandwidth (only the first pole) and to
    % the maximum bandwidth (all the other poles)
    p_m_z = length(p_)-length(z_);
    if p_m_z > 0
        z_(end+1:end+p_m_z) = [exp(-2*pi*tgt_bw*Ts) repmat(exp(-2*pi*max_bw*Ts),1,p_m_z-1)];
    end
    
    % Build filter transfer function and guarantee unitary DC gain
    F{i} = zpk(p_,z_,1,Ts);
    F{i} = F{i}/dcgain(F{i});
    
    % Reduce filter order if specified
    nf(i) = order(F{i});
    if n > 0
        F{i} = balred(F{i},min(order(F{i}),n));
    end
    
    % Compute combined response (actuator and filter)
    AF{i} = A{i}*F{i};
end
