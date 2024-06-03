function p = bwlpred(a, dly, Ts)
% BWLPRED Band-limited discrete-time predictor.
%
% p = BWLPRED(a, dly, Ts)
%
% INPUTS:
%   a:   Parameter to tune predictor bandwidth.
%   dly: Delay (in samples unit) to be predicted (multiple of 0.5).
%   Ts:  (Optional input) Sample time. If not specified -1 is assumed.
%
% OUTPUTS:
%   p:   Transfer function object of the band-limited predictor.

if a < 1
    error('''a'' must be equal or greater than 1 to result in a stable predictor.');
end

N = 2*dly;
if rem(N,1) ~= 0 || N < 0
    error('Delay ''dly'' must be positive and multiple of 0.5');
end

if nargin < 3 || isempty(Ts)
    Ts = -1;
end

p = tf([1+a 0], [a 1], Ts)^N;